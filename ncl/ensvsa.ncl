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
  a = addfile("../sens/ensvsa-TE-m1_a.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  d = 6
  ug       = a->ugrd(d,:,:,:)
  vg       = a->vgrd(d,:,:,:)
  temp     = a->t(d,:,:,:)
  spfh	   = a->q(d,:,:,:)
  ps       = a->pres_meansealev(d,:,:)
  te       = a->te(d,:,:)
  lat_uv   = a->lat
  lon_uv   = a->lon
  plev     = a->lev
  time     = a->time
  
  wks = gsn_open_wks("png","ensvsa_m_t+wind_a_100912_v")
  plot = new(5,graphic)
  cmap = read_colormap_file("posneg_1")
    
  vdate  = cd_string(time(d), "%Y%N%D%H")
  do k = 0,2
     u = ug(k,:,:)
     v = vg(k,:,:)
     t = temp(k,:,:)
     q = spfh(k,:,:)
     q = q*1000
     hgt = plev(k)
;---Create plot
     res                = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = "T,q+wind_"+hgt+"hPa"
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnLeftString	    = "K"
     res@gsnRightString   = "10^(-3)kg/kg"    
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False              ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@gsnScalarContour     = True               ; contours desired
     res@cnFillPalette        = cmap

     ;res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpFillOn             = False
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
;     res@cnMaxLevelCount	   = 10
   if (d.eq.0) then
;  d=0
      res@cnMinLevelValF       = -3.               ; set min contour level
      res@cnMaxLevelValF       = 3.                 ; set max contour level
      res@cnLevelSpacingF      = 0.5               ; set contour spacing
   else if (d.eq.6) then 
;  d=6
      res@cnMinLevelValF       = -5.               ; set min contour level
      res@cnMaxLevelValF       = 5.                 ; set max contour level
      res@cnLevelSpacingF      = 1.               ; set contour spacing
   end if
   end if
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
  
      plott = gsn_csm_vector_scalar_map_ce(wks,u,v,t,res)

      res@cnFillOn             = False             ; turn on color for contours
      res@cnLinesOn            = True              ; turn off contour lines
      res@cnLineLabelsOn       = True
      res@cnLineColor          = "blue"
      res@gsnContourZeroLineThicknessF = 0. 	  ; doubles thickness of zero contour
      res@gsnContourNegLineDashPattern = 1 	  ; sets negative contours to dash pattern 1
      res@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
      res@cnMinLevelValF       = -1.5              ; set min contour level
      res@cnMaxLevelValF       = 1.5                 ; set max contour level
      res@cnLevelSpacingF      = 0.25               ; set contour spacing
; when using the gsn_csm plot templates and overlay, it is necessary to
; turn off the automatic label strings so that you don't get text on 
; top of text
      res@gsnRightString   = ""    
      res@gsnLeftString    = ""
      res@gsnCenterString  = ""
      
      plotq = gsn_csm_contour(wks,q,res)
      overlay(plott,plotq)
      plot(k) = plott
      delete(plott)
      delete(plotq)
   end do

   resps                = True               ; plot mods desired
   resps@gsnDraw	      = False
   resps@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False
     
   resps@gsnCenterString	   = "slp"
   resps@gsnCenterStringOrthogonalPosF = 0.02
   resps@gsnLeftString	    = "hPa"
   resps@gsnAddCyclic	      = False
   resps@cnFillOn             = True               ; turn on color for contours
   resps@cnLinesOn            = True              ; turn off contour lines
   resps@cnLineLabelsOn       = False              ; turn off contour line labels
   resps@gsnScalarContour     = True               ; contours desired
   resps@cnFillPalette        = cmap

   resps@mpLandFillColor      = "gray"            ; set land to be gray
   resps@mpMinLonF            =  120.               ; select a subregion
   resps@mpMaxLonF            =  160.
   resps@mpMinLatF            =  10. 
   resps@mpMaxLatF            =  60.

   resps@lbOrientation            = "Vertical"     ; vertical label bar
   resps@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

   resps@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
   if (d.eq.0) then
;  d=0
      resps@cnMinLevelValF       = -5.0               ; set min contour level
      resps@cnMaxLevelValF       = 5.                 ; set max contour level
      resps@cnLevelSpacingF      = 0.5               ; set contour spacing
   else if (d.eq.6) then 
;  d=6
      resps@cnMinLevelValF       = -10.0               ; set min contour level
      resps@cnMaxLevelValF       = 10.                 ; set max contour level
      resps@cnLevelSpacingF      = 1.               ; set contour spacing
   end if
   end if
   
   plot(3) = gsn_csm_contour_map(wks, ps, resps)

   delete(cmap)
   cmap = read_colormap_file("spread_15lev")
   reste                = True               ; plot mods desired
   reste@gsnDraw	      = False
   reste@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False
     
   reste@gsnCenterString	   = "total moist energy"
   reste@gsnCenterStringOrthogonalPosF = 0.02
   reste@gsnLeftString	    = "J/kg/m^2"
   reste@gsnAddCyclic	      = False
   reste@cnFillOn             = True               ; turn on color for contours
   reste@cnLinesOn            = True              ; turn off contour lines
   reste@cnLineLabelsOn       = False              ; turn off contour line labels
   reste@gsnScalarContour     = True               ; contours desired
   reste@cnFillPalette        = cmap

   reste@mpLandFillColor      = "gray"            ; set land to be gray
   reste@mpMinLonF            =  120.               ; select a subregion
   reste@mpMaxLonF            =  160.
   reste@mpMinLatF            =  10. 
   reste@mpMaxLatF            =  60.

   reste@lbOrientation            = "Vertical"     ; vertical label bar
   reste@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

   reste@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
   if (d.eq.0) then
;  d=0
      reste@cnMinLevelValF       = 0.0               ; set min contour level
      reste@cnMaxLevelValF       = 5.0                 ; set max contour level
      reste@cnLevelSpacingF      = 0.5               ; set contour spacing
   else if (d.eq.6) then 
;  d=6
      reste@cnMinLevelValF       = 0.0               ; set min contour level
      reste@cnMaxLevelValF       = 20.0                 ; set max contour level
      reste@cnLevelSpacingF      = 2.               ; set contour spacing
   end if
   end if
   plot(4) = gsn_csm_contour_map(wks, te, reste)

   resP	                   = True
   resP@gsnPanelMainString	   = "ensvsa-TE T&wind,slp,TE "+vdate
   gsn_panel(wks,plot,(/2,3/),resP)
end
