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
  a = addfile("ensvsa-TE-m1.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  ug       = a->ugrd(6,:,:,:)
  vg       = a->vgrd(6,:,:,:)
  temp     = a->t(6,:,:,:)
  ps       = a->pres_meansealev(6,:,:)
  te       = a->te(6,:,:)
  lat_uv   = a->lat
  lon_uv   = a->lon
  plev     = a->lev
  
  wks = gsn_open_wks("png","ensvsa_t+wind_100912_v")
  plot = new(5,graphic)
  cmap = read_colormap_file("posneg_1")
    
  vdate = 2019101212
  do k = 0,2
     u = ug(k,:,:)
     v = vg(k,:,:)
     t = temp(k,:,:)
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

     res@mpLandFillColor      = "gray"            ; set land to be gray
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
      res@cnMinLevelValF       = -1.0               ; set min contour level
      res@cnMaxLevelValF       = 1                 ; set max contour level
      res@cnLevelSpacingF      = 0.2               ; set contour spacing

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
   resps@cnMinLevelValF       = -3.0               ; set min contour level
   resps@cnMaxLevelValF       = 3                 ; set max contour level
   resps@cnLevelSpacingF      = 0.2               ; set contour spacing

   plot(3) = gsn_csm_contour_map(wks, ps, resps)

   delete(cmap)
   cmap = read_colormap_file("spread_15lev")
   reste                = True               ; plot mods desired
   reste@gsnDraw	      = False
   reste@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False
     
   reste@gsnCenterString	   = "total dry energy"
   reste@gsnCenterStringOrthogonalPosF = 0.02
   reste@gsnLeftString	    = ""
   reste@gsnAddCyclic	      = False
   reste@cnFillOn             = True               ; turn on color for contours
   reste@cnLinesOn            = True              ; turn off contour lines
   reste@cnLineLabelsOn       = False              ; turn off contour line labels
   reste@gsnScalarContour     = True               ; contours desired
   reste@cnFillPalette        = cmap

   reste@mpLandFillColor      = "gray"            ; set land to be gray
   reste@mpMinLonF            =  120.               ; select a subregion
   reste@mpMaxLonF            =  150.
   reste@mpMinLatF            =  10. 
   reste@mpMaxLatF            =  60.

   reste@lbOrientation            = "Vertical"     ; vertical label bar
   reste@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

   reste@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
   reste@cnMinLevelValF       = 0.0               ; set min contour level
   reste@cnMaxLevelValF       = 3.0                 ; set max contour level
   reste@cnLevelSpacingF      = 0.2               ; set contour spacing

   plot(4) = gsn_csm_contour_map(wks, te, reste)

   resP	                   = True
   resP@gsnPanelMainString	   = "ensvsa-TE T&wind,slp,TE "+vdate
   gsn_panel(wks,plot,(/2,3/),resP)
end
