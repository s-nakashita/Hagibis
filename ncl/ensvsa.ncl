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
  a = addfile("../sens/ensvsa-TE-m1-"+orig+"-"+yyyymmddhh+"_n.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  d = 0
  ug       = a->ugrd(d,:,:,:)
  vg       = a->vgrd(d,:,:,:)
  temp     = a->t(d,:,:,:)
  spfh	   = a->q(d,:,:,:)
  ps       = a->pres_meansealev(d,:,:)
; 00Z -> *(-1)
  ug = (/-1*ug/)
  vg = (/-1*vg/)
  temp = (/-1*temp/)
  spfh = (/-1*spfh/)
  ps = (/-1*ps/)
  te       = a->te(d,:,:)
  lat_uv   = a->lat
  lon_uv   = a->lon
  plev     = a->lev
  time     = a->time
  
;  out = "ensvsa_00Zv";"ensvsa_m1_t+wind_n_100912"
;  dev = "png"
  mmddhh = str_get_cols(yyyymmddhh, 4, 9)
  wks = gsn_open_wks("png","ensvsa-TE_m1_t+wind_n_"+orig+"_"+mmddhh)
  plot = new(6,graphic)
  cmap = read_colormap_file("hotcold_18lev")
    
  vdate  = cd_string(time(d), "%Y%N%D%H")

;************************************************
; create points for box
;************************************************
 ypts = (/ 37.0, 37.0, 33.0, 33.0, 37.0/)
 xpts = (/137.0,142.0,142.0,137.0,137.0/)
;************************************************
; add the box
;************************************************
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "black"                     ; color of lines
  resp@gsLineThicknessF = 3.0                       ; thickness of lines
;  resp@gsLineLabelString= "test"                    ; adds a line label string

  do k = 0,2
     u = ug(k,:,:)
     v = vg(k,:,:)
     t = temp(k,:,:)
     q = spfh(k,:,:)
     q = q*1000
     hgt = plev(k)
     ;wks = gsn_open_wks(dev,out+"_"+hgt)
;---Create plot
     res                = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = "T,q+wind_"+hgt+"hPa"
     res@gsnCenterStringFontHeightF = 0.025
     res@gsnCenterStringOrthogonalPosF = 0.04
     ;res@gsnCenterStringParallelPosF   = 0.45
     res@gsnLeftString	    = "K"
     res@gsnLeftStringFontHeightF = 0.02
     res@gsnRightString   = "10^(-3)kg/kg"
     res@gsnRightStringFontHeightF = 0.02
     res@gsnRightStringParallelPosF = 1.1    
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
   if (d.eq.0) then
     res@mpMinLatF            =  0. 
     res@mpMaxLatF            =  60.
   else
     res@mpMinLatF            =  20. 
     res@mpMaxLatF            =  60.
   end if
     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

      res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
   if (d.eq.0) then
;  d=0
      res@cnMinLevelValF       = -2.              ; set min contour level
      res@cnMaxLevelValF       = 2.                ; set max contour level
      res@cnLevelSpacingF      = .2               ; set contour spacing
      res@vcRefMagnitudeF      = 10.0
   else if (d.ge.6) then 
;  d=6
      res@cnMinLevelValF       = -20.               ; set min contour level
      res@cnMaxLevelValF       = 20.                 ; set max contour level
      res@cnLevelSpacingF      = 2.               ; set contour spacing
      res@vcRefMagnitudeF      = 100.0
   end if
   end if
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
      if (d.eq.0) then
;  d=0
      res@cnMinLevelValF       = -5.               ; set min contour level
      res@cnMaxLevelValF       = 5.                 ; set max contour level
      res@cnLevelSpacingF      = 1.               ; set contour spacing
   else if (d.ge.6) then 
;  d=6
      res@cnMinLevelValF       = -50.               ; set min contour level
      res@cnMaxLevelValF       = 50.                 ; set max contour level
      res@cnLevelSpacingF      = 5.               ; set contour spacing
   end if
   end if
; when using the gsn_csm plot templates and overlay, it is necessary to
; turn off the automatic label strings so that you don't get text on 
; top of text
      res@gsnRightString   = ""    
      res@gsnLeftString    = ""
      res@gsnCenterString  = ""
      
      plotq = gsn_csm_contour(wks,q,res)
      overlay(plott,plotq)
; create array of dummy graphic variables. This is required, b/c each line
; must be associated with a unique dummy variable.
      dum = new(4,graphic)  
; draw each line separately. Each line must contain two points.
      do i = 0 , 3
         dum(i)=gsn_add_polyline(wks,plott,xpts(i:i+1),ypts(i:i+1),resp)      
      end do   
      plot(k) = plott
      ;draw(plott)
      ;frame(wks)
      delete(plott)
      delete(plotq)
      delete(dum)
   end do

   delete(cmap)
   cmap = read_colormap_file("spread_15lev")
   ;wks = gsn_open_wks(dev,out+"_te")
   reste                = True               ; plot mods desired
   reste@gsnDraw	      = False
   reste@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False
     
   reste@gsnCenterString	   = "Total moist energy + slp"
   reste@gsnCenterStringOrthogonalPosF = 0.04
   ;reste@gsnCenterStringParallelPosF   = 0.55
   reste@gsnCenterStringFontHeightF    = 0.025
   reste@gsnLeftString	    = "J/kg/m^2"
   reste@gsnLeftStringFontHeightF = 0.02
   reste@gsnLeftStringParallelPosF = -0.1
   reste@gsnRightString	    = "hPa"
   reste@gsnRightStringFontHeightF = 0.02
   reste@gsnRightStringParallelPosF = 1.1
   reste@gsnAddCyclic	      = False
   reste@cnFillOn             = True               ; turn on color for contours
   reste@cnLinesOn            = False             ; turn off contour lines
   reste@cnLineLabelsOn       = False              ; turn off contour line labels
   reste@cnFillPalette        = cmap

   reste@mpLandFillColor      = "gray"            ; set land to be gray
   reste@mpMinLonF            =  120.               ; select a subregion
   reste@mpMaxLonF            =  160.
   if (d.eq.0) then
     reste@mpMinLatF            =  0. 
     reste@mpMaxLatF            =  60.
   else
     reste@mpMinLatF            =  20. 
     reste@mpMaxLatF            =  60.
   end if
   reste@lbOrientation            = "Vertical"     ; vertical label bar
   reste@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

   reste@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
   if (d.eq.0) then
;  d=0
      reste@cnMinLevelValF       = 10.               ; set min contour level
      reste@cnMaxLevelValF       = 100.0                 ; set max contour level
      reste@cnLevelSpacingF      = 10.               ; set contour spacing
   else if (d.ge.6) then 
;  d=6
      reste@cnMinLevelValF       = 100.               ; set min contour level
      reste@cnMaxLevelValF       = 1000.0                 ; set max contour level
      reste@cnLevelSpacingF      = 100.               ; set contour spacing
   end if
   end if
   plote = gsn_csm_contour_map(wks, te, reste)

   resps                = True               ; plot mods desired
   resps@gsnDraw	      = False
   resps@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False
     
   resps@gsnCenterString	   = ""
   resps@gsnCenterStringOrthogonalPosF = 0.02
   resps@gsnLeftString	    = ""
   resps@gsnAddCyclic	      = False
   resps@cnFillOn             = False               ; turn on color for contours
   resps@cnLinesOn            = True              ; turn off contour lines
   resps@cnLineLabelsOn       = False              ; turn off contour line labels
   resps@cnLineColor          = "blue"
   resps@cnLineThicknessF     = 2.
   resps@gsnContourZeroLineThicknessF = 0. 	  ; doubles thickness of zero contour
   resps@gsnContourNegLineDashPattern = 1 	  ; sets negative contours to dash pattern 1
      
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
   else if (d.ge.6) then 
;  d=6
      resps@cnMinLevelValF       = -100.0               ; set min contour level
      resps@cnMaxLevelValF       = 100.                 ; set max contour level
      resps@cnLevelSpacingF      = 10.               ; set contour spacing
   end if
   end if
   
   plotp = gsn_csm_contour(wks, ps, resps)
   overlay(plote,plotp)
; create array of dummy graphic variables. This is required, b/c each line
; must be associated with a unique dummy variable.
   dum = new(4,graphic)  
; draw each line separately. Each line must contain two points.
   do i = 0 , 3
      dum(i)=gsn_add_polyline(wks,plote,xpts(i:i+1),ypts(i:i+1),resp)      
   end do   
   plot(4) = plote
   ;draw(plote)
   ;frame(wks)

   resP	                   = True
   resP@gsnPanelMainString	   = "ensvsa-TE T&wind,slp,TE "+vdate+" "+orig
   gsn_panel(wks,plot,(/2,3/),resP)
end
