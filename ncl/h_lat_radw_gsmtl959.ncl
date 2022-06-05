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
  Exp = "clim"
  if (Exp .eq. "clim") then
    fname = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/np_fcst_asia_"+yyyymmddhh+".nc"
  else if (Exp .eq. "est") then
    fname = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/np_fcst_asia_"+yyyymmddhh+".nc"
  else if (Exp .eq. "mgdsst")
    fname = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/np_fcst_asia_"+yyyymmddhh+".nc"
  end if
  end if
  end if
  f     = addfile (fname,"r")
  tw   = f->U(:,:,:,:) ; tangential wind = E-W component in TC-polar coordinate
  rw   = f->V(:,:,:,:) ; radial wind = N-S component in TC-polar coordinate
  omg  = f->OMG(:,:,:,:) ; hPa/h
  p1d  = f->level
  ;print(p1d)
  omg  = omg * 1e2 / 3600. ; hPa/h -> Pa/s
  omg@units = "Pa/s"
  p    = conform(omg, p1d, 1)
  p    = p*100
  p@units = "Pa"
  temp = f->T(:,:,:,:)
  vw   = omega_to_w(omg, p, temp) ; vertical wind
  delete(p)
  delete(omg)
  delete(temp)
  printMinMax(vw, 0)
  t    = f->time
  ntmax = dimsizes(t)
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
  ncol = 4
  plot_num = nrow*ncol
  plot = new(plot_num,graphic)
  ;ft = tofloat(t)
  ;ft@units = units
  ;ft@calendar = t@calendar
  date  = cd_string(t, "%Y%N%D%H")
  date6h  = cd_string(t6h, "%Y%N%D%H")
  print(date6h)
  t_offset = 0
  wks   = gsn_open_wks ("png", "windrv_h-lat_"+Exp+"_"+date(t_offset)+"-"+date(min((/t_offset+plot_num-1,ntmax-1/))) )          ; send graphics to PNG file
  ;wks   = gsn_open_wks ("png", "windrv_h-lat_6havg_"+Exp+"_"+date6h(t_offset)+"-"+date6h(min((/t_offset+plot_num-2,ntmax6h-1/))) )          ; send graphics to PNG file
  r = new(dimsizes(rw&lat), typeof(rw&lat))
  deg2r(rw&lat, r)
  r@units = "km"
  ;print(r)
  print(sprintf("%.1f", r))
do d=0,plot_num-1
  ind_t2 = d
  ;if (d+t_offset .gt. ntmax6h-1) then
  ;  break
  ;end if
  ;ind_t2 = ind(t.eq.t6h(d+t_offset))
  ;if (d+t_offset.eq.0) then
  ;  ind_t1 = 1
  ;else
  ;  ind_t1 = ind(t.eq.t6h(d+t_offset-1))
  ;end if
  u   = tw(ind_t2,:,:,:)
  ;u   = dim_avg_n(tw(ind_t1:ind_t2,:,:,:),0)
  v   = rw(ind_t2,:,:,:)
  ;v   = dim_avg_n(rw(ind_t1:ind_t2,:,:,:),0)
  w   = vw(ind_t2,:,:,:)
  ;w   = dim_avg_n(vw(ind_t1:ind_t2,:,:,:),0)
  time  = date(ind_t2)
  print(time)
  ;date  = cd_calendar(time, -3)
  ;date  = cd_string(time, "%Y%N%D%H")
  u&level@units = "hPa"
  v&level@units = "hPa"
  w&level@units = "hPa"
  u_zm = u(:,:,0)
  ;printVarSummary(u_zm)
  u_zm = zonalAve(u)
  ;printVarSummary(u_zm)
  v_zm = v(:,:,0)
  ;printVarSummary(v_zm)
  v_zm = zonalAve(v)
  ;printVarSummary(v_zm)
  w_zm = w(:,:,0)
  ;printVarSummary(w_zm)
  w_zm = zonalAve(w)
  ;printVarSummary(w_zm)

  res1                      = True                 ; plot mods desired
  res1@gsnDraw	             = False
  res1@gsnFrame	           = False
  res1@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res1@cnLevelSpacingF      =  3.                  ; contour spacing
  res1@cnMinLevelValF       =  -30.                 ; min level
  res1@cnMaxLevelValF       =  30.                 ; max level
  ;res1@cnMonoLineThickness  = False
  ;res1@cnLineThicknesses    = (/1.0,1.0,1.0,1.0,3.0,1.0,1.0,1.0,1.0,3.0,1.0,1.0,1.0,1.0,3.0/)
  ;res1@cnMonoLineDashPattern= False
  ;res1@cnLineDashPatterns   = (/  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
  res1@cnLinesOn            = False
  res1@cnLineLabelsOn       = False                ; turn off line labels
  res1@cnFillOn             = True                 ; turn on color fill
  res1@cnFillPalette        = "BlWhRe"             ; choose colormap
  res1@tiYAxisString        = "pressure (" + u&level@units + ")"   
  res1@tmYROn               = False
  res1@tmYRLabelsOn         = False
  res1@tmYUseLeft           = False
  res1@trYReverse           = True
  ;res@tmXBMode             = "Explicit"
  ;res@tmXBValues           = r
  ;res@tmXBLabels           = r
  res1@tiMainString         = date(ind_t2)
  
  res2 = True
  res2@gsnDraw = False
  res2@gsnFrame = False
  res2@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res2@cnLevelSpacingF      =  .1                  ; contour spacing
  res2@cnMinLevelValF       =  -1.5                ; min level
  res2@cnMaxLevelValF       =  1.5                 ; max level
  res2@gsnContourNegLineDashPattern = 2
  res2@gsnContourZeroLineThicknessF = 0.0
; Note: u is already on pressure levels. If this were model data, it
; would be necessary to interpolate from the hybrid coordinates to 
; pressure levels.
  data = v_zm * -1.0 ; inflow(northward)=negative, outflow(southward)=positive
  data!0 = "level"
  data&level = u&level
  data!1 = "distance"
  data&distance = r
  data@standard_name = "Radial wind"
  res1@gsnLeftString       = "Radial wind (shade)"
  res1@gsnRightString      = ""
  ;plot(d)  = gsn_csm_pres_hgt(wks, u_zm({100:1000},{75:90}), res )   ; place holder
  plot(d)  = gsn_csm_contour(wks, data({100:1000},{0:500}), res1 )   ; place holder
  data2 = w_zm 
  data2!0 = "level"
  data2&level = w&level
  data2!1 = "distance"
  data2&distance = r
  ;data2@standard_name = "Tangential wind"
  data2@standard_name = "Vertical wind"
  res2@gsnLeftString = ""
  res2@gsnRightString = "Vertical wind (contour)"
  plot2 = gsn_csm_contour(wks, data2({100:1000},{0:500}), res2 )
  overlay(plot(d), plot2)
  delete(data)
  delete(data2)
end do  
  gsn_panel(wks, plot, (/nrow,ncol/), False)
end
