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
;  pres00   = a->HGT(2,:,:,:)
;  pres12   = b->HGT(0,:,:,:)
  u00      = a->UGRD(:,0:2,:,:)
  v00      = a->VGRD(:,0:2,:,:)
  u12      = b->UGRD(:,0:2,:,:)
  v12      = b->VGRD(:,0:2,:,:)
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

  ds = 6
  de = 11
  vdate = 2019100912
  year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
  month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
  day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
  hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
  sec = cd_inv_calendar(year,month,day,hour,0,0, time_12@units, 0)
  newsec = sec
  newsec = sec + 6*60*60*ds
  vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
  ;do l = 0,4 ; level
  wks = gsn_open_wks("png","wind_0912-0900_vint_"+date_12(ds)+"to"+date_12(de));+level(l)+"hPa")
  plot = new(6,graphic)
  cmap = read_colormap_file("sunshine_9lev")
    
  do d = 0,5
     print(vdate)
     ind_00   = ind(date_00.eq.vdate)
     ind_12   = ind(date_12.eq.vdate)
     print(ind_00)
     print(ind_12)

     u = u00(ind_00,:,:,:)
     v = v00(ind_00,:,:,:)
     ;pres = pres12(ind_12,:,:)
     ;pres = pres12(ind_12,:,:) - pres00(ind_00,:,:)
     ndim = dimsizes(u)
     do ilev=0,ndim(0)-1
      do ilat=0,ndim(1)-1
        do ilon=0,ndim(2)-1
          u(ilev,ilat,ilon) = u(ilev,ilat,ilon)*sig(ilev)
          v(ilev,ilat,ilon) = v(ilev,ilat,ilon)*sig(ilev)
        end do
      end do
     end do
     up00 = dim_sum_n(u,0)
     up00!0 = "lat"
     up00!1 = "lon"
     up00&lat = lat_uv
     up00&lon = lon_uv
     up00@_FillValue = default_fillvalue("float")
     vp00 = up00
     vp00 = dim_sum_n(v,0)
     delete(u)
     delete(v)
     u = u12(ind_12,:,:,:)
     v = v12(ind_12,:,:,:)
     ndim = dimsizes(u)
     do ilev=0,ndim(0)-1
      do ilat=0,ndim(1)-1
        do ilon=0,ndim(2)-1
          u(ilev,ilat,ilon) = u(ilev,ilat,ilon)*sig(ilev)
          v(ilev,ilat,ilon) = v(ilev,ilat,ilon)*sig(ilev)
        end do
      end do
     end do
     up12 = dim_sum_n(u,0)
     up12!0 = "lat"
     up12!1 = "lon"
     up12&lat = lat_uv
     up12&lon = lon_uv
     up12@_FillValue = default_fillvalue("float")
     vp12 = up12
     vp12 = dim_sum_n(v,0)
     delete(u)
     delete(v)
     u = up12
     u = up12 - up00
     v = vp12
     v = vp12 - vp00
     w = u
     w = sqrt(u^2+v^2)
     printVarSummary(w)
     ;printMinMax(pres,1)
;---Create plot
     res                      = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = vdate
     res@gsnLeftString	   = ""
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = True               ; turn off contour lines
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
     res@cnMaxLevelCount	   = 10
     res@cnMinLevelValF       = 2.0               ; set min contour level
     res@cnMaxLevelValF       = 20.0                 ; set max contour level
     res@cnLevelSpacingF      = 2.0               ; set contour spacing

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

     u = where(w.lt.1, u@_FillValue, u)
  
     plot(d) = gsn_csm_vector_scalar_map_ce(wks,u,v,w,res)

     year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
     month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
     day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
     hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
     sec = cd_inv_calendar(year,month,day,hour,0,0, time_12@units, 0)
     newsec = sec
     newsec = sec + 6*60*60
     vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
     delete(u)
     delete(v)
   end do
   resP	                   = True
   resP@gsnPanelMainString	   = "wind difference 100912-100900 vertical integrated";@"+level(l)+"hPa"
   gsn_panel(wks,plot,(/2,3/),resP)
   delete(plot)
   ;end do
end
