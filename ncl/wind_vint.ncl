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
  init   = 2019100900
  yyyymmddhh = sprinti("%0.10d",init)
  yyyy = str_get_cols(sprinti("%0.10d",init), 0, 3)
  mmddhh = str_get_cols(sprinti("%0.10d",init), 4, 9)
  ;a = addfile("../../netcdf/tigge/2019/jma/"+mmddhh+"_mean.nc","r")
  fu = addfile("../../netcdf/tigge/2019/jma/rotate/inv_glb_"+yyyymmddhh+"_mean_ug.nc","r")
  fv = addfile("../../netcdf/tigge/2019/jma/rotate/inv_glb_"+yyyymmddhh+"_mean_vg.nc","r")
;  b = addfile("../../netcdf/tigge/2019/jma/100912_mean.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
;  pres00   = a->HGT(2,:,:,:)
;  z00    = a->HGT(:,:,:,:)
  u00    = fu->var(:,2:,:,:) ; 300, 500, 850hPa
  v00    = fv->var(:,2:,:,:)
  lat_uv = fu->lat
  lon_uv = fu->lon
  level  = fu->level
  time   = fu->time
  date   = cd_string(time, "%Y%N%D%H")
  print(date)
  
  sig = new(3,float)
  sig(2) = 3./7. + 4./49.
  sig(1) = 10./49. + 1./7.
  sig(0) = 1./7.
  print(sig)
  
  ds = ind(date.eq.yyyymmddhh)
  de = ind(date.eq."2019101212")
  nd = de - ds + 1
  vdate = init
  year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
  month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
  day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
  hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
  sec = cd_inv_calendar(year,month,day,hour,0,0, time@units, 0)
  newsec = sec
  newsec = sec + 6;*60*60*ds
  vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
  wks = gsn_open_wks("png","asymwind_vint_"+mmddhh+"_"+date(ds)+"to"+date(de))
  plot = new(nd,graphic)
  cmap = read_colormap_file("sunshine_9lev")
    
  do d = 0,nd-1
     print(vdate)
     id   = ind(date.eq.vdate)
     print(id)

     u = u00(id,:,:,:)
     v = v00(id,:,:,:)
     u = u * conform_dims(dimsizes(u), sig, 0)
     v = v * conform_dims(dimsizes(v), sig, 0)
     ndim = dimsizes(u)
     ;do ilev=0,ndim(0)-1
     ; do ilat=0,ndim(1)-1
     ;   do ilon=0,ndim(2)-1
     ;     u(ilev,ilat,ilon) = u(ilev,ilat,ilon)*sig(ilev)
     ;     v(ilev,ilat,ilon) = v(ilev,ilat,ilon)*sig(ilev)
     ;   end do
     ; end do
     ;end do
     up = dim_sum_n(u,0)
     up!0 = "lat"
     up!1 = "lon"
     up&lat = lat_uv
     up&lon = lon_uv
     up@_FillValue = default_fillvalue("float")
     vp = up
     vp = dim_sum_n(v,0)
     w = up
     w = sqrt(up^2+vp^2)
     printVarSummary(w)
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
     res@mpMinLatF            =  -10. 
     res@mpMaxLatF            =  30.

     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.
     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     ;res@cnMaxLevelCount	     = 10
     res@cnMinLevelValF       = 5.0               ; set min contour level
     res@cnMaxLevelValF       = 50.0                 ; set max contour level
     res@cnLevelSpacingF      = 5.0               ; set contour spacing

     res@vcRefMagnitudeF           = 10.0             ; define vector ref mag
     res@vcRefLengthF              = 0.045           ; define length of vec ref
     res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
     res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
     res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

;     res@vcMinMagnitudeF	      = 1.0		  ; no effect with CurlyVector
     res@vcGlyphStyle            = "CurlyVector"     ; turn on curly vectors
     res@vcLineArrowColor        = "black"           ; change vector color
     res@vcLineArrowThicknessF   = 1.0               ; change vector thickness
     res@vcMinDistanceF	        = 0.02
     res@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last

     ;up = where(w.lt.5., up@_FillValue, up)
  
     plot(d) = gsn_csm_vector_scalar_map(wks,up,vp,w,res)
     
     year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
     month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
     day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
     hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
     sec = cd_inv_calendar(year,month,day,hour,0,0, time@units, 0)
     newsec = sec
     newsec = sec + 6;*60*60
     vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
     delete(u)
     delete(v)
     delete(res)
   end do
   resP	                   = True
   resP@gsnPanelMainString	   = "asymmetric vertical integrated wind "+mmddhh
   gsn_panel(wks,plot,(/4,(nd+1)/4/),resP)
   delete(plot)
end
