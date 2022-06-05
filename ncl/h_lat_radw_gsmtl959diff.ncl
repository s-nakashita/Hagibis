;----------------------------------------------------------------------
; h_lat_6.ncl
;----------------------------------------------------------------------
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
;----------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "calc_e.ncl"
begin
  ;f     = addfile ("../../netcdf/gsm/gl/2019/10/rotate/np_init_vort.nc","r")
  yyyymmddhh = 2019100912
  f1 = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/np_fcst_asia_"+yyyymmddhh+".nc","r")
  ;Exp = "clim"
  ;if (Exp .eq. "clim") then
  f2 = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/np_fcst_asia_"+yyyymmddhh+".nc","r")
  ;else if (Exp .eq. "mgdsst")
  ;  f2 = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/np_fcst_asia_"+yyyymmddhh+".nc","r")
  f3 = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/np_fcst_asia_"+yyyymmddhh+".nc","r")
  ;end if
  ;end if
  t    = f1->time
  ntmax = dimsizes(t)
  tw1   = f1->U(:,:,:,:) ; tangential wind = E-W component in TC-polar coordinate
  rw1   = f1->V(:,:,:,:) ; radial wind = N-S component in TC-polar coordinate
  tw2   = f2->U(:,:,:,:) ; tangential wind = E-W component in TC-polar coordinate
  rw2   = f2->V(:,:,:,:) ; radial wind = N-S component in TC-polar coordinate
  tw3   = f3->U(:,:,:,:) ; tangential wind = E-W component in TC-polar coordinate
  rw3   = f3->V(:,:,:,:) ; radial wind = N-S component in TC-polar coordinate
  p1d  = f1->level
  ;print(p1d)
  p    = conform(tw1, p1d, 1)
  p    = p*100
  p@units = "Pa"
  omg  = f1->OMG(:,:,:,:)
  omg  = omg * 1e2 / 3600.
  omg@units = "Pa/s"
  temp = f1->T(:,:,:,:)
  vw1   = omega_to_w(omg, p, temp) ; vertical wind
  delete(omg)
  delete(temp)
  omg  = f2->OMG(:,:,:,:)
  omg  = omg * 1e2 / 3600.
  omg@units = "Pa/s"
  temp = f2->T(:,:,:,:)
  vw2   = omega_to_w(omg, p, temp) ; vertical wind
  delete(omg)
  delete(temp)
  omg  = f3->OMG(:,:,:,:)
  omg  = omg * 1e2 / 3600.
  omg@units = "Pa/s"
  temp = f3->T(:,:,:,:)
  vw3   = omega_to_w(omg, p, temp) ; vertical wind
  delete(omg)
  delete(temp)
  delete(p)
  f6h = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/np_fcst_phy_"+yyyymmddhh+".nc","r")
  t6h = f6h->time
  ntmax6h = dimsizes(t6h)
  delete(f6h)
  ;yyyy = 2019
  ;mm   = 10
  ;dd   = 12
  ;hh   = 12
  ;units= t@units
  ;opt  = 0
  ;opt@calendar = t@calendar
  ;vt   = cd_inv_calendar(yyyy, mm, dd, hh, 0, 0, units, opt)
  ;print(vt)
  ;id   = ind(t.eq.vt)
  ;print(id)
  ;plot_num = id/2 + 1
  nrow = 3
  ncol = 3
  plot_num = nrow*ncol
  plot = new(plot_num,graphic)
  ;ft = tofloat(t)
  ;ft@units = units
  ;ft@calendar = t@calendar
  date    = cd_string(t, "%Y%N%D%H")
  date6h  = cd_string(t6h, "%Y%N%D%H")
  print(date6h)
  t_offset = 4
  ;wks   = gsn_open_wks ("png", "winds_h-lat_"+Exp+"-est_"+date(t_offset)+"-"+date(min((/t_offset+plot_num-1,ntmax-1/))) )          ; send graphics to PNG file
  ;wks   = gsn_open_wks ("png", "winds_h-lat_6havg_"+Exp+"-mgdsst_"+date6h(t_offset)+"-"+date6h(min((/t_offset+plot_num-2,ntmax6h-1/))) )          ; send graphics to PNG file
  wks   = gsn_open_wks ("png", "winds_h-lat_diff_"+date(t_offset)+"-"+date(min((/t_offset+nrow-1,ntmax-1/))) )          ; send graphics to PNG file
  r = new(dimsizes(rw1&lat), typeof(rw1&lat))
  deg2r(rw1&lat, r)
  r@units = "km"
  print(r)
  print(sprintf("%.1f", r))
do d=0,nrow-1
  ;if (d+t_offset .gt. ntmax6h-1) then
  ;  break
  ;end if
  ;ind_t2 = ind(t.eq.t6h(d+t_offset))
  ;if (d+t_offset.eq.0) then
  ;  ind_t1 = 1
  ;else
  ;  ind_t1 = ind(t.eq.t6h(d+t_offset-1))
  ;end if
  ind_t2 = d+t_offset
  u1   = tw1(ind_t2,:,:,:)
  ;u1   = dim_avg_n(tw1(ind_t1:ind_t2,:,:,:),0)
  v1   = rw1(ind_t2,:,:,:)
  ;v1   = dim_avg_n(rw1(ind_t1:ind_t2,:,:,:),0)
  w1   = vw1(ind_t2,:,:,:)
  ;w1   = dim_avg_n(vw1(ind_t1:ind_t2,:,:,:),0)
  u2   = tw2(ind_t2,:,:,:)
  ;u2   = dim_avg_n(tw2(ind_t1:ind_t2,:,:,:),0)
  v2   = rw2(ind_t2,:,:,:)
  ;v2   = dim_avg_n(rw2(ind_t1:ind_t2,:,:,:),0)
  w2   = vw2(ind_t2,:,:,:)
  ;w2   = dim_avg_n(vw2(ind_t1:ind_t2,:,:,:),0)
  u3   = tw3(ind_t2,:,:,:)
  ;u2   = dim_avg_n(tw2(ind_t1:ind_t2,:,:,:),0)
  v3   = rw3(ind_t2,:,:,:)
  ;v2   = dim_avg_n(rw2(ind_t1:ind_t2,:,:,:),0)
  w3   = vw3(ind_t2,:,:,:)
  ;w2   = dim_avg_n(vw2(ind_t1:ind_t2,:,:,:),0)
  ;u = u1
  u2 = u2 - u1
  ;u = u1 - u2
  ;v = v1
  v2 = v2 - v1
  ;v = v1 - v2
  ;w = w1
  w2 = w2 - w1
  ;w = w1 - w2
  u3 = u3 - u1
  v3 = v3 - v1
  w3 = w3 - w1
  time  = t(ind_t2)
  print(time)
  ;date  = cd_calendar(time, -3)
  ;date  = cd_string(time, "%Y%N%D%H")
  ;u&level@units = "hPa"
  ;v&level@units = "hPa"
  ;w&level@units = "hPa"
  ;u_zm = u(:,:,0)
  ;printVarSummary(u_zm)
  ;u_zm = zonalAve(u)
  ;printVarSummary(u_zm)
  ;v_zm = v(:,:,0)
  ;printVarSummary(v_zm)
  ;v_zm = zonalAve(v)
  ;printVarSummary(v_zm)
  ;w_zm = w(:,:,0)
  ;printVarSummary(w_zm)
  ;w_zm = zonalAve(w)
  ;printVarSummary(w_zm)
  u1_zm = u1(:,:,0)
  u1_zm = zonalAve(u1)
  v1_zm = v1(:,:,0)
  v1_zm = zonalAve(v1)
  w1_zm = w1(:,:,0)
  w1_zm = zonalAve(w1)
  u2_zm = u2(:,:,0)
  u2_zm = zonalAve(u2)
  v2_zm = v2(:,:,0)
  v2_zm = zonalAve(v2)
  w2_zm = w2(:,:,0)
  w2_zm = zonalAve(w2)
  u3_zm = u3(:,:,0)
  u3_zm = zonalAve(u3)
  v3_zm = v3(:,:,0)
  v3_zm = zonalAve(v3)
  w3_zm = w3(:,:,0)
  w3_zm = zonalAve(w3)

  res1                      = True                 ; plot mods desired
  res1@gsnDraw	             = False
  res1@gsnFrame	           = False
  res1@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  ;res1@cnLevelSpacingF      =  1.                  ; contour spacing
  ;res1@cnMinLevelValF       =  -10.                 ; min level
  ;res1@cnMaxLevelValF       =  10.                 ; max level
  ;res1@cnMonoLineThickness  = False
  ;res1@cnLineThicknesses    = (/1.0,1.0,1.0,1.0,3.0,1.0,1.0,1.0,1.0,3.0,1.0,1.0,1.0,1.0,3.0/)
  ;res1@cnMonoLineDashPattern= False
  ;res1@cnLineDashPatterns   = (/  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
  res1@cnLinesOn            = False
  res1@cnLineLabelsOn       = False                ; turn off line labels
  res1@cnFillOn             = True                 ; turn on color fill
  res1@cnFillPalette        = "BlWhRe"             ; choose colormap
  res1@tiYAxisString        = "pressure (" + u1&level@units + ")"   
  res1@tmYROn               = False
  res1@tmYRLabelsOn         = False
  res1@tmYUseLeft           = False
  res1@trYReverse           = True
  ;res@tmXBMode             = "Explicit"
  ;res@tmXBValues           = r
  ;res@tmXBLabels           = r
  res1@gsnRightString       = u1@units
  ;res1@tiMainString         = date(ind_t2)
  
  res2 = True
  res2@gsnDraw = False
  res2@gsnFrame = False
  res2@gsnRightString = ""
  res2@gsnLeftString = ""
  res2@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  ;res2@cnLevelSpacingF      =  .05                  ; contour spacing
  ;res2@cnMinLevelValF       =  -0.5                 ; min level
  ;res2@cnMaxLevelValF       =  .5                 ; max level
  res2@gsnContourNegLineDashPattern = 2
  res2@gsnContourZeroLineThicknessF = 0.0
; Note: u is already on pressure levels. If this were model data, it
; would be necessary to interpolate from the hybrid coordinates to 
; pressure levels.
  ; EST : row data
  data = v1_zm * -1.0 ; inflow(northward)=negative, outflow(southward)=positive
  data!0 = "level"
  data&level = u1&level
  data!1 = "distance"
  data&distance = r
  data@standard_name = "Radial wind"
  res1@gsnLeftString       = "Radial wind (shade)"
  res1@gsnRightString      = ""
  res1@cnLevelSpacingF     =  3.                  ; contour spacing
  res1@cnMinLevelValF      =  -30.                 ; min level
  res1@cnMaxLevelValF      =  30.                 ; max level
  res1@tiMainString        = date(ind_t2) + ", EST"
  ;plot(d)  = gsn_csm_pres_hgt(wks, u_zm({100:1000},{75:90}), res )   ; place holder
  plot(3*d)  = gsn_csm_contour(wks, data({100:1000},{0:500}), res1 )   ; place holder
  data2 = w1_zm 
  data2!0 = "level"
  data2&level = u1&level
  data2!1 = "distance"
  data2&distance = r
  ;data2@standard_name = "Tangential wind"
  data2@standard_name = "Vertical wind"
  res2@gsnLeftString = ""
  res2@gsnRightString = data2@standard_name+" (contour)"
  res2@cnLevelSpacingF      =  .1                  ; contour spacing
  res2@cnMinLevelValF       =  -1.5                 ; min level
  res2@cnMaxLevelValF       =  1.5                 ; max level
  plot2 = gsn_csm_contour(wks, data2({100:1000},{0:500}), res2 )
; Add shading
  opt = True
  opt@gsnShadeFillType = "pattern"
  opt@gsnShadeHigh     = 17
  plot2 = gsn_contour_shade(plot2, -0.1, 10.1, opt)
  delete(opt)
  overlay(plot(3*d), plot2)
  delete(data)
  delete(data2)
  ; CLIM : difference
  data = v2_zm * -1.0 ; inflow(northward)=negative, outflow(southward)=positive
  data!0 = "level"
  data&level = u1&level
  data!1 = "distance"
  data&distance = r
  data@standard_name = "Radial wind"
  res1@gsnLeftString       = "Radial wind (shade)"
  res1@gsnRightString      = ""
  res1@cnLevelSpacingF     =  1.                  ; contour spacing
  res1@cnMinLevelValF      =  -10.                 ; min level
  res1@cnMaxLevelValF      =  10.                 ; max level
  res1@tiMainString        = "CLIM - EST"
  ;plot(d)  = gsn_csm_pres_hgt(wks, u_zm({100:1000},{75:90}), res )   ; place holder
  plot(3*d+1)  = gsn_csm_contour(wks, data({100:1000},{0:500}), res1 )   ; place holder
  data2 = w2_zm 
  data2!0 = "level"
  data2&level = u1&level
  data2!1 = "distance"
  data2&distance = r
  ;data2@standard_name = "Tangential wind"
  data2@standard_name = "Vertical wind"
  res2@gsnLeftString = ""
  res2@gsnRightString = data2@standard_name+" (contour)"
  res2@cnLevelSpacingF      =  .05                  ; contour spacing
  res2@cnMinLevelValF       =  -.5                 ; min level
  res2@cnMaxLevelValF       =  .5                 ; max level
  plot2 = gsn_csm_contour(wks, data2({100:1000},{0:500}), res2 )
; Add shading
  opt = True
  opt@gsnShadeFillType = "pattern"
  opt@gsnShadeHigh     = 17
  plot2 = gsn_contour_shade(plot2, -0.1, 10.1, opt)
  delete(opt)
  overlay(plot(3*d+1), plot2)
  delete(data)
  delete(data2)
  ; MGD : difference
  data = v3_zm * -1.0 ; inflow(northward)=negative, outflow(southward)=positive
  data!0 = "level"
  data&level = u1&level
  data!1 = "distance"
  data&distance = r
  data@standard_name = "Radial wind"
  res1@gsnLeftString       = "Radial wind (shade)"
  res1@gsnRightString      = ""
  res1@cnLevelSpacingF     =  1.                  ; contour spacing
  res1@cnMinLevelValF      =  -10.                 ; min level
  res1@cnMaxLevelValF      =  10.                 ; max level
  res1@tiMainString        = "MGD - EST"
  ;plot(d)  = gsn_csm_pres_hgt(wks, u_zm({100:1000},{75:90}), res )   ; place holder
  plot(3*d+2)  = gsn_csm_contour(wks, data({100:1000},{0:500}), res1 )   ; place holder
  data2 = w3_zm 
  data2!0 = "level"
  data2&level = u1&level
  data2!1 = "distance"
  data2&distance = r
  ;data2@standard_name = "Tangential wind"
  data2@standard_name = "Vertical wind"
  res2@gsnLeftString = ""
  res2@gsnRightString = data2@standard_name+" (contour)"
  res2@cnLevelSpacingF      =  .05                  ; contour spacing
  res2@cnMinLevelValF       =  -.5                 ; min level
  res2@cnMaxLevelValF       =  .5                 ; max level
  plot2 = gsn_csm_contour(wks, data2({100:1000},{0:500}), res2 )
; Add shading
  opt = True
  opt@gsnShadeFillType = "pattern"
  opt@gsnShadeHigh     = 17
  plot2 = gsn_contour_shade(plot2, -0.1, 10.1, opt)
  delete(opt)
  overlay(plot(3*d+2), plot2)
  delete(data)
  delete(data2)
end do  
  gsn_panel(wks, plot, (/nrow,ncol/), False)
end
