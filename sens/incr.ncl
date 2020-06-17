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
  a = addfile("../../netcdf/tigge/2019/jma/100900_ctrl.nc","r")
  b = addfile("../../netcdf/tigge/2019/jma/100912_ctrl.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  u00      = a->UGRD(:,:,:,:)
  v00      = a->VGRD(:,:,:,:)
  t00      = a->TMP(:,:,:,:)
  ps00     = a->PRES_meansealevel(:,:,:)
  u12      = b->UGRD(:,:,:,:)
  v12      = b->VGRD(:,:,:,:)
  t12      = b->TMP(:,:,:,:)
  ps12     = b->PRES_meansealevel(:,:,:)
  lat_uv   = a->lat
  lon_uv   = a->lon
  plev     = a->level
  time_00   = a->time
  time_12   = b->time
  date_00  = cd_string(time_00, "%Y%N%D%H")
  date_12  = cd_string(time_12, "%Y%N%D%H")
  
  wks = gsn_open_wks("png","incr_100912-100900_c")
  plot = new(4,graphic)
  cmap = read_colormap_file("hotcold_18lev")
    
  vdate = 2019100912
  ind_00 = ind(date_00.eq.vdate)
  ind_12 = ind(date_12.eq.vdate)
  do k = 0,2
     u = u12(ind_12,k,:,:)
     u = u12(ind_12,k,:,:) - u00(ind_00,k,:,:) 
     v = v12(ind_12,k,:,:)
     v = v12(ind_12,k,:,:) - v00(ind_00,k,:,:)
     t = t12(ind_12,k,:,:)
     t = t12(ind_12,k,:,:) - t00(ind_00,k,:,:)
     ps = ps12(ind_12,:,:)
     ps = ps12(ind_12,:,:) - ps00(ind_00,:,:)
     hgt = plev(k)
;---Create plot
     res                = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = "T+wind_"+hgt+"hPa"
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnLeftString	    = ""
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False              ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@gsnScalarContour     = True               ; contours desired
     res@cnFillPalette        = cmap

;     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpFillOn	      = False
     res@mpMinLonF            =  120.               ; select a subregion
     res@mpMaxLonF            =  150.
     res@mpMinLatF            =  10. 
     res@mpMaxLatF            =  60.

     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

      res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
      res@cnMinLevelValF       = -3.0               ; set min contour level
      res@cnMaxLevelValF       = 3                 ; set max contour level
      res@cnLevelSpacingF      = 0.5               ; set contour spacing

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
      u = where(wind.lt.0.1, u@_FillValue, u)
  
      plot(k) = gsn_csm_vector_scalar_map_ce(wks,u,v,t,res)
   end do

   resps                = True               ; plot mods desired
   resps@gsnDraw	      = False
   resps@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False
     
   resps@gsnCenterString	   = "slp"
   resps@gsnCenterStringOrthogonalPosF = 0.02
   resps@gsnLeftString	    = ""
   resps@gsnAddCyclic	      = False
   resps@cnFillOn             = True               ; turn on color for contours
   resps@cnLinesOn            = True              ; turn off contour lines
   resps@cnLineLabelsOn       = False              ; turn off contour line labels
   resps@gsnScalarContour     = True               ; contours desired
   resps@cnFillPalette        = cmap

   resps@mpLandFillColor      = "gray"            ; set land to be gray
   resps@mpMinLonF            =  120.               ; select a subregion
   resps@mpMaxLonF            =  150.
   resps@mpMinLatF            =  10. 
   resps@mpMaxLatF            =  60.

   resps@lbOrientation            = "Vertical"     ; vertical label bar
   resps@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

   resps@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
   resps@cnMinLevelValF       = -500.0               ; set min contour level
   resps@cnMaxLevelValF       = 500.0                 ; set max contour level
   resps@cnLevelSpacingF      = 50.0               ; set contour spacing

   plot(3) = gsn_csm_contour_map(wks, ps, resps)

   resP	                   = True
   resP@gsnPanelMainString	   = "increment(ctrl) T&wind,slp "+vdate
   gsn_panel(wks,plot,(/2,2/),resP)
end
