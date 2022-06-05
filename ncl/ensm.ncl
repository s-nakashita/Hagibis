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
  orig = "jma"
  dev = "png"
  init = "0900"
  if (init.eq."0900") then
    d_offset = 2
  else
    d_offset = 0
  end if
  a = addfile("/Volumes/dandelion/netcdf/tigge/2019/"+orig+"/glb_201910"+init+"_mean.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  um   = a->u(:,:,:,:)
  vm   = a->v(:,:,:,:)
  tm   = a->t(:,:,:,:)
  qm   = a->q(:,:,:,:)
  psm  = a->msl(:,:,:)
  psm  = psm*0.01
  zm   = a->gh(:,:,:,:)
  lat  = a->latitude
  lon  = a->longitude
  plev = a->level
  time = a->time
  date = cd_string(time, "%Y%N%D%H")
  lon0 = 120.
  lon1 = 180.
  lat0 = -10.
  lat1 = 60.
  
  out = "ensm_"+orig+"_"+init
  wks = gsn_open_wks(dev,out+"_slp")
  plot = new(4,graphic)
  cmap = read_colormap_file("CBR_coldhot")
  
  resps                = True               ; plot mods desired
  resps@gsnDraw	      = False
  resps@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False
  resps@gsnCenterStringOrthogonalPosF = 0.02
  resps@gsnLeftString	    = ""
  resps@gsnRightString	    = "hPa"
  resps@gsnRightStringFontHeightF = 0.02
  resps@gsnAddCyclic	      = False

  resps@mpLandFillColor      = "gray"            ; set land to be gray
  resps@mpMinLonF            =  lon0               ; select a subregion
  resps@mpMaxLonF            =  lon1
  resps@mpMinLatF            =  lat0
  resps@mpMaxLatF            =  lat1

  resps@lbOrientation            = "Vertical"     ; vertical label bar
  ;resps@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

  resps@vpXF = 0.15
  resps@vpWidthF = 0.7

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

  resps@cnFillOn             = False             ; turn on color for contours
  resps@cnLinesOn            = True              ; turn off contour lines
  resps@cnLineLabelsOn       = True
  resps@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
   ;resps@cnMaxLevelCount	   = 10
  resps@cnMinLevelValF       = 980.0              ; set min contour level
  resps@cnMaxLevelValF       = 1020.0                 ; set max contour level
  resps@cnLevelSpacingF      = 4.0               ; set contour spacing

  i = 0
  do d = 0, 12, 4
  resps@gsnCenterString	   = date(d+d_offset)
  plot(i) = gsn_csm_contour_map(wks,psm(d+d_offset,:,:),resps)
  i = i + 1
  end do
  gsn_panel(wks,plot,(/2,2/),False)
  delete(plot)
  delete(wks)
  do k = 3, 8
    u = um(:,k,:,:)
    v = vm(:,k,:,:)
    t = tm(:,k,:,:)
    q = qm(:,k,:,:)
    ;q = q*10^3
    z = zm(:,k,:,:)
    hgt = plev(k)
    ;;; calculate equivalent potential temperature
    ;; spfh -> rh
    rh = relhum(t, q, conform(t, hgt*100.0, -1))
    ;; t, rh -> tlcl
    tlcl = tlcl_rh_bolton(t,rh,(/1,1/))
    ;; theta_e
    theta_e = pot_temp_equiv_tlcl(conform(t, hgt*100.0, -1), t, tlcl, q, (/1,1,0,1/))
    printMinMax(theta_e,0)
    plot = new(4,graphic)
    wks = gsn_open_wks(dev,out+"_"+hgt)
;---Create plot
    res                = True               ; plot mods desired
    res@gsnDraw	      = False
    res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
    res@gsnCenterStringOrthogonalPosF = 0.02
    res@gsnLeftString	    = "gpm"
    res@gsnLeftStringFontHeightF = 0.02
    res@gsnRightString	    = "K"
    res@gsnRightStringFontHeightF = 0.02
    res@gsnRightStringParallelPosF = 1.1    
    res@gsnAddCyclic	      = False
    res@cnFillOn             = True               ; turn on color for contours
    res@cnLinesOn            = False              ; turn off contour lines
    res@cnLineLabelsOn       = False              ; turn off contour line labels
    res@gsnScalarContour     = True               ; contours desired
    res@cnFillPalette        = cmap

;    res@mpLandFillColor      = "gray"            ; set land to be gray
    res@mpFillOn	      = False
    res@mpMinLonF            =  lon0               ; select a subregion
    res@mpMaxLonF            =  lon1
    res@mpMinLatF            =  lat0
    res@mpMaxLatF            =  lat1

    res@lbOrientation            = "Vertical"     ; vertical label bar
    res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer
    res@vpXF = 0.15
    res@vpWidthF = 0.7

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

    res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;;     res@cnMaxLevelCount	   = 10
    if(k.eq.3) then ;300hPa
      res@cnMinLevelValF       = 320.0             ; set min contour level
      res@cnMaxLevelValF       = 350.0             ; set max contour level
      res@cnLevelSpacingF      = 2.                ; set contour spacing
      res@vcRefMagnitudeF      = 40.0              ; define vector ref mag
    else if (k.eq.4) then ;500hPa
      res@cnMinLevelValF       = 302.0             ; set min contour level
      res@cnMaxLevelValF       = 332.0             ; set max contour level
      res@cnLevelSpacingF      = 2.                ; set contour spacing
      res@vcRefMagnitudeF      = 30.0              ; define vector ref mag
    else if (k.eq.5) then ;700hPa
      res@cnMinLevelValF       = 292.0             ; set min contour level
      res@cnMaxLevelValF       = 322.0             ; set max contour level
      res@cnLevelSpacingF      = 2.                ; set contour spacing
      res@vcRefMagnitudeF      = 25.0              ; define vector ref mag
    else if (k.eq.6) then ;850hPa
      res@cnMinLevelValF       = 282.0             ; set min contour level
      res@cnMaxLevelValF       = 312.0             ; set max contour level
      res@cnLevelSpacingF      = 2.                ; set contour spacing
      res@vcRefMagnitudeF      = 25.0              ; define vector ref mag
    else ;925,1000hPa
      res@cnMinLevelValF       = 276.0             ; set min contour level
      res@cnMaxLevelValF       = 304.0             ; set max contour level
      res@cnLevelSpacingF      = 2.                ; set contour spacing
      res@vcRefMagnitudeF      = 25.0              ; define vector ref mag
    end if
    end if
    end if
    end if
    res@vcRefLengthF              = 0.045           ; define length of vec ref
    res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector

;     res@vcMinMagnitudeF	      = 1.0		  ; no effect with CurlyVector
    ;res@vcGlyphStyle            = "CurlyVector"     ; turn on curly vectors
    res@vcFillArrowsOn       = True
    res@vcMinFracLengthF     = 0.2
    res@vcFillArrowEdgeColor = "black"
    res@vcFillArrowFillColor = "green"
    res@vcMinDistanceF	      = 0.02
    res@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last

    wind = sqrt(u^2+v^2)
    u = where(wind.lt.0.5, u@_FillValue, u)
  
    ;plott = gsn_csm_vector_scalar_map_ce(wks,u,v,t,res)
    ;plott = gsn_csm_contour_map_ce(wks,t,res)

    resz = True
    resz@gsnDraw	      = False
    resz@gsnFrame	      = False
    resz@cnFillOn             = False             ; turn on color for contours
    resz@cnLinesOn            = True              ; turn off contour lines
    resz@cnLineLabelsOn       = True
    resz@cnLineColor          = "black"
    resz@cnLineThicknessF     = 3.0
    ;resz@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
    ;;res@cnMaxLevelCount	   = 10
    ;if (k .eq. 5) then
    ;   res@cnMinLevelValF       = 1200.0              ; set min contour level
    ;   res@cnMaxLevelValF       = 1500.0                 ; set max contour level
    ;   res@cnLevelSpacingF      = 30.0               ; set contour spacing
    ;end if
    ;if (k .eq. 4) then
    ;   res@cnMinLevelValF       = 5200.0              ; set min contour level
    ;   res@cnMaxLevelValF       = 5800.0                 ; set max contour level
    ;   res@cnLevelSpacingF      = 50.0               ; set contour spacing
    ;end if
    ;if (k .eq. 3) then
    ;   res@cnMinLevelValF       = 8800.0              ; set min contour level
    ;   res@cnMaxLevelValF       = 9600.0                 ; set max contour level
    ;   res@cnLevelSpacingF      = 80.0               ; set contour spacing
    ;end if
; when using the gsn_csm plot templates and overlay, it is necessary to
; turn off the automatic label strings so that you don't get text on 
; top of text
    resz@gsnRightString   = ""    
    resz@gsnLeftString    = ""
    resz@gsnCenterString  = ""

    i = 0
    do d=0,12,4
      res@gsnCenterString = date(d+d_offset)
      plot(i) = gsn_csm_vector_scalar_map_ce(wks,u(d+d_offset,:,:),v(d+d_offset,:,:),\
      theta_e(d+d_offset,:,:),res)
      printMinMax(z(d+d_offset,{lat0:lat1},{lon0:lon1}),0)
      plotz = gsn_csm_contour(wks,z(d+d_offset,:,:),resz)
      overlay(plot(i),plotz)
      delete(plotz)
      i = i + 1
    end do
    gsn_panel(wks,plot,(/2,2/),False)
    delete(plot)
    delete(wks)
   end do

;   resP	                   = True
;   resP@gsnPanelMainString	   = "perturbation member 26 T&Q&wind,slp "+vdate
;   gsn_panel(wks,plot,(/1,4/),resP)
end
