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
  yyyymmddhh = 2019100900
  Exp = "clim"
  if (Exp .eq. "clim") then
    fname  = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/np_fcst_asia_"+yyyymmddhh+".nc"
    frname = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/np_fcst_phy_"+yyyymmddhh+".nc"
  else if (Exp .eq. "est") then
    fname  = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/np_fcst_asia_"+yyyymmddhh+".nc"
    frname = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/np_fcst_phy_"+yyyymmddhh+".nc"
  else if (Exp .eq. "mgdsst")
    fname  = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/np_fcst_asia_"+yyyymmddhh+".nc"
    frname = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/np_fcst_phy_"+yyyymmddhh+".nc"
  end if
  end if
  end if
  f    = addfile (fname,"r")
  fr   = addfile (frname,"r")
  tw   = f->U(:,:,:,:) ; tangential wind = E-W component in TC-polar coordinate
  hrcv = fr->HRCV(:,:,:,:)
  hrlc = fr->HRLC(:,:,:,:)
  ;hrr  = hrcv
  hrr  = hrlc
  ;hrr  = hrcv + hrlc
  printMinMax(hrr, 0)
  t    = f->time
  tr   = fr->time
  ntmax = dimsizes(tr)
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
  date  = cd_string(tr, "%Y%N%D%H")
  print(date)
  t_offset = 0
  wks   = gsn_open_wks ("png", "tanw+hrlc_h-lat_"+Exp+"_"+date(t_offset)+"-"+date(min((/t_offset+plot_num-1,ntmax-1/))) )          ; send graphics to PNG file
  r = new(dimsizes(tw&lat), typeof(tw&lat))
  deg2r(tw&lat, r)
  r@units = "km"
  print(sprintf("%.1f", r))
  r2 = new(dimsizes(hrr&lat), typeof(hrr&lat))
  deg2r(hrr&lat, r2)
  r2@units = "km"
  print(sprintf("%.1f", r2))
ip = 0
do d=0,ntmax-1
  if (ip .gt. plot_num-1) then
    break
  end if
  ind_t2 = ind(t.eq.tr(d+t_offset))
  if(d+t_offset.eq.0) then
    ind_t1 = 0
  else
    ind_t1 = ind(t.eq.tr(d+t_offset-1))
  end if
  u    = tw(ind_t2,:,:,:)
  u    = dim_avg_n(tw(ind_t1:ind_t2,:,:,:), 0)
  hr   = hrr(d+t_offset,:,:,:)
  time = tr(d+t_offset)
  print(time)
  ;date  = cd_calendar(time, -3)
  ;date  = cd_string(time, "%Y%N%D%H")
  u&level@units = "hPa"
  hr&level@units = "hPa"
  u_zm = u(:,:,0)
  u_zm = zonalAve(u)
  printVarSummary(u_zm)
  hr_zm = hr(:,:,0)
  hr_zm = zonalAve(hr)
  printVarSummary(hr_zm)

  res1                      = True                 ; plot mods desired
  res1@gsnDraw	             = False
  res1@gsnFrame	           = False
  res1@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res1@cnLevelSpacingF      =  1.                  ; contour spacing
  res1@cnMinLevelValF       =  -9.                 ; min level
  res1@cnMaxLevelValF       =  9.                 ; max level
  ;res1@cnMonoLineThickness  = False
  ;res1@cnLineThicknesses    = (/1.0,1.0,1.0,1.0,3.0,1.0,1.0,1.0,1.0,3.0,1.0,1.0,1.0,1.0,3.0/)
  ;res1@cnMonoLineDashPattern= False
  ;res1@cnLineDashPatterns   = (/  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
  res1@cnLinesOn            = False
  res1@cnLineLabelsOn       = False                ; turn off line labels
  res1@cnFillOn             = True                 ; turn on color fill
  res1@cnFillPalette        = "BlueDarkRed18"             ; choose colormap
  res1@tiYAxisString        = "pressure (" + u&level@units + ")"   
  res1@tmYROn               = False
  res1@tmYRLabelsOn         = False
  res1@tmYUseLeft           = False
  res1@trYReverse           = True
  ;res@tmXBMode             = "Explicit"
  ;res@tmXBValues           = r
  ;res@tmXBLabels           = r
  res1@tiMainString         = date(d+t_offset)
  
  res2 = True
  res2@gsnDraw = False
  res2@gsnFrame = False
  res2@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res2@cnLevelSpacingF      =  10.                  ; contour spacing
  res2@cnMinLevelValF       =  5.                 ; min level
  res2@cnMaxLevelValF       =  45.                 ; max level
; Note: u is already on pressure levels. If this were model data, it
; would be necessary to interpolate from the hybrid coordinates to 
; pressure levels.
  data = hr_zm
  data!0 = "level"
  data&level = u&level
  data!1 = "distance"
  data&distance = r2
  data@standard_name = "Heating rate (cloud)"
  res1@gsnLeftString       = "Heating rate (K/h)"
  res1@gsnRightString      = ""
  ;plot(d)  = gsn_csm_pres_hgt(wks, u_zm({100:1000},{75:90}), res )   ; place holder
  plot(ip)  = gsn_csm_contour(wks, data({100:1000},{0:500}), res1 )   ; place holder
  data2 = u_zm 
  data2!0 = "level"
  data2&level = u&level
  data2!1 = "distance"
  data2&distance = r
  data2@standard_name = "Tangential wind"
  res2@gsnLeftString = ""
  res2@gsnRightString = "Tangential wind (contour)"
  plot2 = gsn_csm_contour(wks, data2({100:1000},{0:500}), res2 )
  overlay(plot(ip), plot2)
  delete(data)
  delete(data2)
  ip = ip + 1
end do  
  gsn_panel(wks, plot, (/nrow,ncol/), False)
end
