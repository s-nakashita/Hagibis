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
  z00      = a->HGT(:,:,:,:)
  z12      = b->HGT(:,:,:,:)
  lat_uv   = b->lat
  lon_uv   = b->lon
  level    = b->level
  time_00  = a->time
  time_12  = b->time
  date_00  = cd_string(time_00, "%Y%N%D%H")
  date_12  = cd_string(time_12, "%Y%N%D%H")
  print(date_00)
  print(date_12)

  sig = new(3,float)
  sig(0) = 3./7. + 4./49.
  sig(1) = 10./49. + 1./7.
  sig(2) = 1./7.
  print(sig)

  do l = 0,4 ; level
  ds = 0
  de = 5
  nd = de - ds + 1
  vdate = 2019100912
  year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
  month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
  day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
  hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
  sec = cd_inv_calendar(year,month,day,hour,0,0, time_12@units, 0)
  newsec = sec
  newsec = sec + 6*60*60*ds
  vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
  wks = gsn_open_wks("png","hgt_0912-0900_"+level(l)+"hPa_"+date_12(ds)+"to"+date_12(de))
  plot = new(nd,graphic)
  cmap = read_colormap_file("BlueDarkRed18")
    
  do d = 0,nd-1
     print(vdate)
     ind_00   = ind(date_00.eq.vdate)
     ind_12   = ind(date_12.eq.vdate)
     print(ind_00)
     print(ind_12)

     dz = z00(ind_00,l,:,:)
     dz = z12(ind_12,l,:,:) - z00(ind_00,l,:,:)
     dz@_FillValue = default_fillvalue("float")
     printVarSummary(dz)
     ;printMinMax(pres,1)
;---Create plot
     res                      = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = sprinti("%0.10d",vdate)
     res@gsnLeftString	   = ""
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False               ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@gsnScalarContour     = True               ; (for gsn_csm_vector_scalar_map)draw contours
     res@cnFillPalette        = cmap

     res@mpFillOn             = False 
;     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpMinLonF            =  120.               ; select a subregion
     res@mpMaxLonF            =  160.
     res@mpMinLatF            =  10. 
     res@mpMaxLatF            =  60.

     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMaxLevelCount	   = 20
     res@cnMinLevelValF       = -50.0               ; set min contour level
     res@cnMaxLevelValF       = 50.0                 ; set max contour level
     res@cnLevelSpacingF      = 5.0               ; set contour spacing

     plotdz = gsn_csm_contour_map_ce(wks,dz,res)
     z = z12(ind_12,l,:,:)
     ;if (l.eq.0) then
     ; res@cnMinLevelValF       = 1200.0               ; set min contour level
     ; res@cnMaxLevelValF       = 1500.0                 ; set max contour level
     ; res@cnLevelSpacingF      = 30.0               ; set contour spacing
     ;else if (l.eq.1) then 
     ; res@cnMinLevelValF       = 5200.0               ; set min contour level
     ; res@cnMaxLevelValF       = 5800.0                 ; set max contour level
     ; res@cnLevelSpacingF      = 50.0               ; set contour spacing
     ;else if (l.eq.2) then
     ; res@cnMinLevelValF       = 8800.0               ; set min contour level
     ; res@cnMaxLevelValF       = 9600.0                 ; set max contour level
     ; res@cnLevelSpacingF      = 80.0               ; set contour spacing
     ;else if (l.eq.3) then
     ; res@cnMinLevelValF       = 9600.0               ; set min contour level
     ; res@cnMaxLevelValF       = 11000.0                 ; set max contour level
     ; res@cnLevelSpacingF      = 100.0               ; set contour spacing
     ;else if (l.eq.4) then
     ; res@cnMinLevelValF       = 11000.0               ; set min contour level
     ; res@cnMaxLevelValF       = 13000.0                 ; set max contour level
     ; res@cnLevelSpacingF      = 100.0               ; set contour spacing
     ;end if
     ;end if
     ;end if
     ;end if
     ;end if
     resz                      = True               ; plot mods desired
     resz@gsnDraw	      = False
     resz@gsnFrame	      = False
     resz@cnFillOn             = False               ; turn on color for contours
     resz@cnLinesOn            = True                ; turn off contour lines
     resz@cnLineLabelsOn        = True
     resz@cnLineColor          = "blue"
     resz@gsnRightString   = ""    
     resz@gsnLeftString    = ""
     resz@gsnCenterString  = ""
     plotz = gsn_csm_contour(wks,z,resz)
     overlay(plotdz,plotz)
     plot(d) = plotdz
     delete(plotdz)
     delete(plotz)

     year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
     month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
     day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
     hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
     sec = cd_inv_calendar(year,month,day,hour,0,0, time_12@units, 0)
     newsec = sec
     newsec = sec + 6*60*60
     vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
     delete(dz)
     delete(z)
   end do
   resP	                   = True
   resP@gsnPanelMainString	   = "hgt difference 100912-100900 @"+level(l)+"hPa"
   gsn_panel(wks,plot,(/2,(nd+1)/2/),resP)
   delete(plot)
   end do
end
