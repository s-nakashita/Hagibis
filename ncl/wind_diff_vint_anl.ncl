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
  hh = "12"
  a = addfile("../../netcdf/gsm/gl/2019/10/init.nc","r")
  b = addfile("../../netcdf/tigge/2019/jma/1009"+hh+"_mean.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
;  pres00   = a->HGT(2,:,:,:)
;  pres12   = b->HGT(0,:,:,:)
  ua      = a->UGRD(:,:,180:340,180:360)
  va      = a->VGRD(:,:,180:340,180:360)
  uf      = b->UGRD(:,0:2,:,:)
  vf      = b->VGRD(:,0:2,:,:)
  lat_uv   = b->lat
  lon_uv   = b->lon
  leva    = a->level
  levf    = b->level
  time_a  = a->time
  time_f  = b->time
  date_a  = cd_string(time_a, "%Y%N%D%H")
  date_f  = cd_string(time_f, "%Y%N%D%H")
  print(date_a)
  print(date_f)

  sig = new(3,float)
  sig(0) = 3./7. + 4./49.
  sig(1) = 10./49. + 1./7.
  sig(2) = 1./7.
  print(sig)

  ds = 0
  de = 10
  vdate = 2019100900 + toint(hh)
  year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
  month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
  day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
  hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
  sec = cd_inv_calendar(year,month,day,hour,0,0, time_f@units, 0)
  newsec = sec
  newsec = sec + 6*60*60*ds
  vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
  ;do l = 0,4 ; level
  wks = gsn_open_wks("png","wind_09"+hh+"-anl_vint_"+date_f(ds)+"to"+date_f(de));+level(l)+"hPa")
  plot = new(6,graphic)
  cmap = read_colormap_file("sunshine_9lev")
    
  do d = 0,5
     print(vdate)
     ind_a   = ind(date_a.eq.vdate)
     ind_f   = ind(date_f.eq.vdate)
     print(ind_a)
     print(ind_f)

     u = uf(ind_f,:,:,:)
     v = vf(ind_f,:,:,:)
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
     upf = dim_sum_n(u,0)
     upf!0 = "lat"
     upf!1 = "lon"
     upf&lat = lat_uv
     upf&lon = lon_uv
     upf@_FillValue = default_fillvalue("float")
     vpf = upf
     vpf = dim_sum_n(v,0)
     ;delete(u)
     ;delete(v)
     do l = 0, 2
      il = ind(leva.eq.levf(l))
      u(l,:,:) = (/ua(ind_a,il,:,:)/)
      v(l,:,:) = (/va(ind_a,il,:,:)/)
     end do
     ndim = dimsizes(u)
     do ilev=0,ndim(0)-1
      do ilat=0,ndim(1)-1
        do ilon=0,ndim(2)-1
          u(ilev,ilat,ilon) = u(ilev,ilat,ilon)*sig(ilev)
          v(ilev,ilat,ilon) = v(ilev,ilat,ilon)*sig(ilev)
        end do
      end do
     end do
     upa = dim_sum_n(u,0)
     upa!0 = "lat"
     upa!1 = "lon"
     upa&lat = lat_uv
     upa&lon = lon_uv
     upa@_FillValue = default_fillvalue("float")
     vpa = upa
     vpa = dim_sum_n(v,0)
     delete(u)
     delete(v)
     u = upf
     u = upf - upa
     v = vpf
     v = vpf - vpa
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
     sec = cd_inv_calendar(year,month,day,hour,0,0, time_f@units, 0)
     newsec = sec
     newsec = sec + 12*60*60
     vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
     delete(u)
     delete(v)
   end do
   resP	                   = True
   resP@gsnPanelMainString	   = "wind difference 1009"+hh+"-analysis vertical integrated";@"+level(l)+"hPa"
   gsn_panel(wks,plot,(/2,3/),resP)
   delete(plot)
   ;end do
end
