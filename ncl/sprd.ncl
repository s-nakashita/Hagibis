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
;  orig = "jma"
;  init = "2019100900"
;  regrid = False
;  dev = "png"
  if (regrid) then
  a = addfile("/Volumes/dandelion/netcdf/tigge/2019/regrid/"+orig+"/glb_"+init+"_mean.nc","r")
  b = addfile("/Volumes/dandelion/netcdf/tigge/2019/regrid/"+orig+"/glb_"+init+"_spr.nc","r")
  else
  a = addfile("/Volumes/dandelion/netcdf/tigge/2019/"+orig+"/glb_"+init+"_mean.nc","r")
  b = addfile("/Volumes/dandelion/netcdf/tigge/2019/"+orig+"/glb_"+init+"_spr.nc","r")
  end if
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  ;u00      = a->UGRD(:,:,:,:)
  ;v00      = a->VGRD(:,:,:,:)
  ;t00      = a->TMP(:,:,:,:)
  ;q00      = a->SPFH(:,:,:,:)
  ;ps00     = a->PRES_meansealevel(:,:,:)
  ;z00      = a->HGT(:,:,:,:)
  ;u12      = b->UGRD(:,:,:,:)
  ;v12      = b->VGRD(:,:,:,:)
  ;t12      = b->TMP(:,:,:,:)
  ;q12      = b->SPFH(:,:,:,:)
  ;ps12     = b->PRES_meansealevel(:,:,:)
  ;z12      = b->HGT(:,:,:,:)
  u00      = a->u(:,:,:,:)
  v00      = a->v(:,:,:,:)
  t00      = a->t(:,:,:,:)
  ps00     = a->msl(:,:,:)
  z00      = a->gh(:,:,:,:)
  u12      = b->u(:,:,:,:)
  v12      = b->v(:,:,:,:)
  t12      = b->t(:,:,:,:)
  ps12     = b->msl(:,:,:)
  z12      = b->gh(:,:,:,:)
  if (orig .ne. "ukmo") then
  q00      = a->q(:,:,:,:)
  q12      = b->q(:,:,:,:)
  end if
  if(regrid) then
  lat_uv   = a->lat
  lon_uv   = a->lon
  else
  lat_uv   = a->latitude
  lon_uv   = a->longitude
  end if
  plev     = a->level
  time_00   = a->time
  time_12   = b->time
  date_00  = cd_string(time_00, "%Y%N%D%H")
  date_12  = cd_string(time_12, "%Y%N%D%H")
  lon0 = 110.
  lon1 = 170.
  lat0 = -10.
  lat1 = 60.
  
  if(regrid) then
  out = "sprd_regrid_"+orig+"_"+init
  else
  out = "sprd_"+orig+"_"+init
  end if
;  wks = gsn_open_wks("png","prtb_mem26")
;  plot = new(4,graphic)
  cmap = read_colormap_file("spread_15lev")
    
  vdate = init
  ind_00 = ind(date_00.eq.vdate)
  ind_12 = ind(date_12.eq.vdate)
  print(ind_00)
  print(ind_12)
  ps = ps12(ind_12,:,:)
  ;ps = ps12(ind_12,:,:) - ps00(ind_00,:,:)
  ps = ps * 0.01 ; Pa->hPa
  pful = ps00(ind_12,:,:)
  pful = pful * 0.01 ; Pa->hPa
  wks = gsn_open_wks(dev,out+"_slp")
  resps                = True               ; plot mods desired
  resps@gsnDraw	      = False
  resps@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False
     
  resps@gsnCenterString	   = "slp"
  resps@gsnCenterStringOrthogonalPosF = 0.02
  resps@gsnLeftString	    = ""
  resps@gsnRightString	    = "hPa"
  resps@gsnRightStringFontHeightF = 0.02
  resps@gsnAddCyclic	      = False
  resps@cnFillOn             = True               ; turn on color for contours
  resps@cnLinesOn            = False              ; turn off contour lines
  resps@cnLineLabelsOn       = False              ; turn off contour line labels
  ;resps@gsnScalarContour     = True               ; contours desired
  resps@cnFillPalette        = cmap

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

  resps@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
  resps@cnMinLevelValF       = 0.4              ; set min contour level
  resps@cnMaxLevelValF       = 1.2                 ; set max contour level
  resps@cnLevelSpacingF      = 0.1             ; set contour spacing

  plotp = gsn_csm_contour_map(wks, ps, resps)

  resps@cnFillOn             = False             ; turn on color for contours
  resps@cnLinesOn            = True              ; turn off contour lines
  resps@cnLineLabelsOn       = True
  resps@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
   ;resps@cnMaxLevelCount	   = 10
  resps@cnMinLevelValF       = 980.0              ; set min contour level
  resps@cnMaxLevelValF       = 1028.0                 ; set max contour level
  resps@cnLevelSpacingF      = 4.0               ; set contour spacing
; when using the gsn_csm plot templates and overlay, it is necessary to
; turn off the automatic label strings so that you don't get text on 
; top of text
  resps@gsnRightString   = ""    
  resps@gsnLeftString    = ""
  resps@gsnCenterString  = ""
      
  plotpf = gsn_csm_contour(wks,pful,resps)
  overlay(plotp,plotpf)
  ;plot(0) = plotp
  draw(plotp)
  frame(wks)
  delete(plotp)
  delete(plotpf)
  levels = (/850,500,300/)
  ;do k = 3, 8 ;; jma
  do k = 0, dimsizes(levels)-1 ;; others
    level = levels(k)
     u = u12(ind_12,{level},:,:)
     ;u = u12(ind_12,k,:,:) - u00(ind_00,k,:,:) 
     v = v12(ind_12,{level},:,:)
     ;v = v12(ind_12,k,:,:) - v00(ind_00,k,:,:)
     ke = u 
     ke = 0.5*sqrt(u^2 + v^2)
     ke@long_name = "kinetic energy spread"
     ke@units = "m~S~2~N~/s~S~2~N~"
     t = t12(ind_12,{level},:,:)
     ;t = t12(ind_12,k,:,:) - t00(ind_00,k,:,:)
     if (orig .ne. "ukmo") then
     q = q12(ind_12,{level},:,:)
     ;q = q12(ind_12,k,:,:) - q00(ind_00,k,:,:)
     q = q*10^3
     q@units = "g/ks"
     end if
     z = z12(ind_12,{level},:,:)
     ;z = z12(ind_12,k,:,:) - z00(ind_00,k,:,:)
     wks = gsn_open_wks(dev,out+"_"+level)
     plots = new(2,graphic)
;---Create plot
     res                = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
     res@gsnAddCyclic	      = False
     resc = res

;     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpFillOn	      = False
     res@mpMinLonF            =  lon0               ; select a subregion
     res@mpMaxLonF            =  lon1
     res@mpMinLatF            =  lat0
     res@mpMaxLatF            =  lat1

     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

     res@vpXF = 0.15
     res@vpWidthF = 0.7

;     res@cnInfoLabelOn	      = False
     res1 = res
     res2 = res
     ;res1 : T+q
     res1@gsnCenterString	   = "T,q";+"_"+level+"hPa"
     res1@gsnCenterStringOrthogonalPosF = 0.02
     res1@gsnLeftString	    = "g/kg"
     res1@gsnLeftStringFontHeightF = 0.02
     res1@gsnRightString	    = "K"
     res1@gsnRightStringFontHeightF = 0.02
     res1@gsnRightStringParallelPosF = 1.1    
     res1@cnFillOn             = True               ; turn on color for contours
     res1@cnLinesOn            = False              ; turn off contour lines
     res1@cnLineLabelsOn       = False              ; turn off contour line labels
     ;res@gsnScalarContour     = True               ; contours desired
     res1@cnFillPalette        = cmap

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

      res1@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res1@cnMaxLevelCount	   = 10
      res1@cnMinLevelValF       = 0.4               ; set min contour level
      res1@cnMaxLevelValF       = 1.2                 ; set max contour level
      res1@cnLevelSpacingF      = 0.1               ; set contour spacing

      ;res@vcRefMagnitudeF           = 5.0             ; define vector ref mag
      ;res@vcRefLengthF              = 0.045           ; define length of vec ref
      ;res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
      ;res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
      ;res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

;     res@vcMinMagnitudeF	      = 1.0		  ; no effect with CurlyVector
      ;res@vcGlyphStyle            = "CurlyVector"     ; turn on curly vectors
      ;res@vcLineArrowColor        = "black"           ; change vector color
      ;res@vcLineArrowThicknessF   = 2.0               ; change vector thickness
      ;res@vcMinDistanceF	      = 0.02
      ;res@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last

      ;wind = sqrt(u^2+v^2)
      ;u = where(wind.lt.0.5, u@_FillValue, u)
  
      ;plott = gsn_csm_vector_scalar_map_ce(wks,u,v,t,res)
      plots(0) = gsn_csm_contour_map_ce(wks,t,res1)

     if (orig .ne. "ukmo") then
      resc@cnFillOn             = False             ; turn on color for contours
      resc@cnLinesOn            = True              ; turn off contour lines
      resc@cnLineLabelsOn       = True
      resc@cnLineColor          = "blue"
      resc@cnLineThicknessF     = 3.0
      resc@gsnContourZeroLineThicknessF = 0. 	  ; doubles thickness of zero contour
      resc@gsnContourNegLineDashPattern = 1 	  ; sets negative contours to dash pattern 1
      resc@cnMaxLevelCount	   = 10
;      resc@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
;      resc@cnMinLevelValF       = 0.5             ; set min contour level
;      resc@cnMaxLevelValF       = 3.0                 ; set max contour level
;      resc@cnLevelSpacingF      = 0.5               ; set contour spacing
; when using the gsn_csm plot templates and overlay, it is necessary to
; turn off the automatic label strings so that you don't get text on 
; top of text
      resc@gsnRightString   = ""    
      resc@gsnLeftString    = ""
      resc@gsnCenterString  = ""
      
      plotq = gsn_csm_contour(wks,q,resc)
      overlay(plots(0),plotq)
     end if

      res2@cnFillOn             = True             ; turn on color for contours
      res2@cnLinesOn            = False              ; turn off contour lines
      res2@cnLineLabelsOn       = False
      res2@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
      res2@cnFillPalette        = cmap
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
      res2@cnMinLevelValF = 0.4
      res2@cnMaxLevelValF = 1.2
      res2@cnLevelSpacingF = 0.1
; when using the gsn_csm plot templates and overlay, it is necessary to
; turn off the automatic label strings so that you don't get text on 
; top of text
      res2@gsnLeftString   = ""    
      res2@gsnRightString    = "m~S~2~N~/s~S~2~N~"
      res2@gsnCenterString  = "KE"

      plots(1) = gsn_csm_contour_map_ce(wks,ke,res2)
      
      resP = True
      resP@gsnPanelMainString = level+" hPa"
      gsn_panel(wks,plots,(/1,2/),resP)
      delete(resP)
      delete(res1)
      delete(res2)
      delete(resc)
      delete(res)
   end do

;   resP	                   = True
;   resP@gsnPanelMainString	   = "perturbation member 26 T&Q&wind,slp "+vdate
;   gsn_panel(wks,plot,(/1,4/),resP)
end
