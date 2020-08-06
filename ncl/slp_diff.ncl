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
  ;hh = 12
  a = addfile("../../netcdf/tigge/2019/jma/100900_mean.nc","r")
  b = addfile("../../netcdf/tigge/2019/jma/100912_mean.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  ps00    = a->PRES_meansealevel(:,:,:);PRMSL_meansealevel(:,180:340,180:360)
  ps12    = b->PRES_meansealevel(:,:,:)
  lat_uv  = b->lat
  lon_uv  = b->lon
  time00  = a->time
  time12  = b->time
  date00    = cd_string(time00, "%Y%N%D%H")
  date12    = cd_string(time12, "%Y%N%D%H")
  print(date00)
 
  ds = 6
  de = 11
  nd = de - ds + 1
  vdate = 2019100912 ;+ toint(hh)
  year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
  month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
  day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
  hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
  sec = cd_inv_calendar(year,month,day,hour,0,0, time12@units, 0)
  newsec = sec
  newsec = sec + 6*60*60*ds
  vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
  wks = gsn_open_wks("png","slp_0912-0900_"+date12(ds)+"to"+date12(de))
  plot = new(nd,graphic)
  cmap = read_colormap_file("BlueDarkRed18")
    
  do d = 0,nd-1
     print(vdate)
     ind00   = ind(date00.eq.vdate)
     ind12   = ind(date12.eq.vdate)
     print(ind00)
     print(ind12)
     dz = ps12(ind12,:,:)
     dz = dz - ps00(ind00,:,:)
     
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
     res@cnMinLevelValF       = -500.0               ; set min contour level
     res@cnMaxLevelValF       = 500.0                 ; set max contour level
     res@cnLevelSpacingF      = 50.0               ; set contour spacing

     plotdz = gsn_csm_contour_map_ce(wks,dz,res)
     z = ps12(ind12,:,:)
     resp = True
     resp@gsnDraw	      = False
     resp@gsnFrame	      = False
     ;res@cnMinLevelValF       = 96000.0               ; set min contour level
     ;res@cnMaxLevelValF       = 108000.0                 ; set max contour level
     ;res@cnLevelSpacingF      = 400.0               ; set contour spacing
     resp@cnFillOn             = False               ; turn on color for contours
     resp@cnLinesOn            = True                ; turn off contour lines
     resp@cnLineLabelsOn        = True
     resp@cnLineColor          = "blue"
     resp@gsnRightString   = ""    
     resp@gsnLeftString    = ""
     resp@gsnCenterString  = ""
     plotz = gsn_csm_contour(wks,z,resp)
     overlay(plotdz,plotz)
     plot(d) = plotdz
     delete(plotdz)
     delete(plotz)

     year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
     month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
     day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
     hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
     sec = cd_inv_calendar(year,month,day,hour,0,0, time12@units, 0)
     newsec = sec
     newsec = sec + 6*60*60
     vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
     delete(dz)
     delete(z)
   end do
   resP	                   = True
   resP@gsnPanelMainString	   = "slp difference 100912-100900"
   gsn_panel(wks,plot,(/2,(nd+1)/2/),resP)
   delete(plot)
end
