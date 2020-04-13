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
  a = addfile("../netcdf/jma/2019100900_mean.nc","r")
  b = addfile("../netcdf/jma/2019100912_mean.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  pres00   = a->HGT_200mb(:,:,:)
  pres12   = b->HGT_200mb(:,:,:)
  u00      = a->UGRD_200mb(:,:,:)
  v00      = a->VGRD_200mb(:,:,:)
  u12      = b->UGRD_200mb(:,:,:)
  v12      = b->VGRD_200mb(:,:,:)
  lat_uv   = b->latitude
  lon_uv   = b->longitude
  time_00  = a->time
  time_12  = b->time
  date_00  = cd_string(time_00, "%Y%N%D%H")
  date_12  = cd_string(time_12, "%Y%N%D%H")
  print(date_00)
  print(date_12)

  wks = gsn_open_wks("png","200hgt+wind_100912-100900")
  plot = new(7,graphic)
  cmap = read_colormap_file("ViBlGrWhYeOrRe")
    
  vdate = 2019100912
  do d = 0,6
     print(vdate)
     ind_00   = ind(date_00.eq.vdate)
     ind_12   = ind(date_12.eq.vdate)
     print(ind_00)
     print(ind_12)

     u = u12(ind_12,:,:)
     u = u12(ind_12,:,:) - u00(ind_00,:,:)
     v = v12(ind_12,:,:)
     v = v12(ind_12,:,:) - v00(ind_00,:,:)
     pres = pres12(ind_12,:,:)
     pres = pres12(ind_12,:,:) - pres00(ind_00,:,:)
     printVarSummary(u)
     printMinMax(u,1)
     printVarSummary(v)
     printMinMax(v,1)
     printVarSummary(pres)
     printMinMax(pres,1)
;---Create plot
     res                      = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = "ft="+vdate
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False              ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@gsnScalarContour     = True               ; contours desired
     res@cnFillPalette        = cmap

     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpMinLonF            =  120.               ; select a subregion
     res@mpMaxLonF            =  150.
     res@mpMinLatF            =  0. 
     res@mpMaxLatF            =  60.

     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

      res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
      res@cnMinLevelValF       = -50.0               ; set min contour level
      res@cnMaxLevelValF       = 50                 ; set max contour level
      res@cnLevelSpacingF      = 10.0               ; set contour spacing

      res@vcRefMagnitudeF           = 10.0             ; define vector ref mag
      res@vcRefLengthF              = 0.045           ; define length of vec ref
      res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
      res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
      res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

;     res@vcMinMagnitudeF	      = 1.0		  ; no effect with CurlyVector
      res@vcGlyphStyle            = "CurlyVector"     ; turn on curly vectors
      res@vcLineArrowColor        = "black"           ; change vector color
      res@vcLineArrowThicknessF   = 2.0               ; change vector thickness
      res@vcMinDistanceF	      = 0.02
      res@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last

      wind = sqrt(u^2+v^2)
      u = where(wind.lt.1, u@_FillValue, u)
  
      plot(d) = gsn_csm_vector_scalar_map_ce(wks,u,v,pres,res)

      year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
      month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
      day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
      hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
      sec = cd_inv_calendar(year,month,day,hour,0,0, time_12@units, 0)
      newsec = sec
      newsec = sec + 12*60*60
      vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
   end do
   resP	                   = True
   resP@gsnPanelMainString	   = "hgt&wind difference 100912-100900"
   gsn_panel(wks,plot,(/2,4/),resP)
end
