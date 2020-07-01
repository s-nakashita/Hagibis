;----------------------------------------------------------------------
; vector_1.ncl
;
; Concepts illustrated:
;   - Overlaying vectors and filled contours on a map
;   - Manually attaching lat/lon coordinate arrays to a variable
;   - Changing the length of the smallest vector as a fraction of the reference vector
;   - Moving the vector reference annotation to the top right of the plot
;   - Drawing curly vectors
;   - Setting the color for vectors
;   - Making the labelbar be vertical
;   - Increasing the thickness of vectors
;----------------------------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

begin

;---Open netCDF file
  a = addfile("../../netcdf/tigge/2019/jma/100900_mean.nc","r")
  b = addfile("../../netcdf/tigge/2019/jma/100912_mean.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  time_00  = a->time
  time_12  = b->time
  date_00  = cd_string(time_00, "%Y%N%D%H")
  date_12  = cd_string(time_12, "%Y%N%D%H")
  print(date_00)
  print(date_12)
  vdate = 2019100912
  print(vdate)
  ind_00   = ind(date_00.eq.vdate)
  ind_12   = ind(date_12.eq.vdate)
  print(ind_00)
  print(ind_12)
  hgt00   = a->HGT(ind_00,:,:,:)
  hgt12   = b->HGT(ind_12,:,:,:)
;  u00      = a->UGRD(:,:,:,:)
;  v00      = a->VGRD(:,:,:,:)
;  u12      = b->UGRD(:,:,:,:)
;  v12      = b->VGRD(:,:,:,:)
  lat_uv   = b->lat
  lon_uv   = b->lon
  lev      = b->level

  wks = gsn_open_wks("png","hgt_100912-100900")
  plot = new(5,graphic)
  cmap = read_colormap_file("ViBlGrWhYeOrRe")

  do d = 0,4
     hgt = hgt12(d,:,:)
     dhgt = hgt12(d,:,:)
     dhgt = hgt12(d,:,:) - hgt00(d,:,:)
;---Create plot
     res                      = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = lev(d)+"hPa"
     res@gsnLeftString      = ""
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False              ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@gsnScalarContour     = True               ; contours desired
     res@cnFillPalette        = cmap

     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpMinLonF            =  120.               ; select a subregion
     res@mpMaxLonF            =  160.
     res@mpMinLatF            =  10. 
     res@mpMaxLatF            =  60.

     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
     res@cnMinLevelValF       = -40.0               ; set min contour level
     res@cnMaxLevelValF       = 40.0                 ; set max contour level
     res@cnLevelSpacingF      = 8.0               ; set contour spacing

     plotdz = gsn_csm_contour_map_ce(wks,dhgt,res)
; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

     res@cnFillOn             = False               ; turn on color for contours
     res@cnLinesOn            = True              ; turn off contour lines
     
     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
    if (d.eq.0) then
      res@cnMinLevelValF       = 1200.0               ; set min contour level
      res@cnMaxLevelValF       = 1500.0                 ; set max contour level
      res@cnLevelSpacingF      = 30.0               ; set contour spacing
    else if (d.eq.1) then 
      res@cnMinLevelValF       = 5200.0               ; set min contour level
      res@cnMaxLevelValF       = 5800.0                 ; set max contour level
      res@cnLevelSpacingF      = 50.0               ; set contour spacing
    else if (d.eq.2) then
      res@cnMinLevelValF       = 8800.0               ; set min contour level
      res@cnMaxLevelValF       = 9600.0                 ; set max contour level
      res@cnLevelSpacingF      = 80.0               ; set contour spacing
    else if (d.eq.3) then
      res@cnMinLevelValF       = 9600.0               ; set min contour level
      res@cnMaxLevelValF       = 11000.0                 ; set max contour level
      res@cnLevelSpacingF      = 100.0               ; set contour spacing
    else if (d.eq.4) then
      res@cnMinLevelValF       = 11000.0               ; set min contour level
      res@cnMaxLevelValF       = 13000.0                 ; set max contour level
      res@cnLevelSpacingF      = 150.0               ; set contour spacing
    end if
    end if
    end if
    end if
    end if
     res@cnLineLabelsOn        = True
     res@gsnRightString   = ""    
     res@gsnLeftString    = ""
     res@gsnCenterString  = ""
    
    plotz = gsn_csm_contour(wks,hgt,res)
    overlay(plotdz,plotz)
    plot(d) = plotdz
    ;plot(d) = gsn_csm_contour_map_overlay(wks,dhgt,hgt,res1,res2)
    delete(plotdz)
    delete(plotz)
   end do
   resP	                   = True
   resP@gsnPanelMainString	   = "hgt & hgt difference 100912-100900"
   gsn_panel(wks,plot,(/2,3/),resP)
end
