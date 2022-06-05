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
  f1  = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/np_fcst_asia_"+yyyymmddhh+".nc","r")
  fr1 = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/np_fcst_phy_"+yyyymmddhh+".nc","r")
  Exp = "mgdsst"
  if (Exp .eq. "clim") then
    f2  = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/np_fcst_asia_"+yyyymmddhh+".nc","r")
    fr2 = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/np_fcst_phy_"+yyyymmddhh+".nc","r")
  else if (Exp .eq. "mgdsst")
    f2  = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/np_fcst_asia_"+yyyymmddhh+".nc","r")
    fr2 = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/np_fcst_phy_"+yyyymmddhh+".nc","r")
  end if
  end if
  tw1  = f1->U(:,:,:,:) ; tangential wind = E-W component in TC-polar coordinate
  hrcv = fr1->HRCV(:,:,:,:)
  hrlc = fr1->HRLC(:,:,:,:)
  hrr1 = hrcv
  hrr1 = hrcv + hrlc
  printMinMax(hrr1, 0)
  tw2  = f2->U(:,:,:,:) ; tangential wind = E-W component in TC-polar coordinate
  hrcv = fr2->HRCV(:,:,:,:)
  hrlc = fr2->HRLC(:,:,:,:)
  hrr2 = hrcv
  hrr2 = hrcv + hrlc
  printMinMax(hrr2, 0)
  t    = f1->time
  tr   = fr1->time
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
  wks   = gsn_open_wks ("png", "tanw+hr_h-lat_est-"+Exp+"_"+date(t_offset)+"-"+date(min((/t_offset+plot_num-1,ntmax-1/))) )          ; send graphics to PNG file
  r = new(dimsizes(tw1&lat), typeof(tw1&lat))
  deg2r(tw1&lat, r)
  r@units = "km"
  print(sprintf("%.1f", r))
  r2 = new(dimsizes(hrr1&lat), typeof(hrr1&lat))
  deg2r(hrr1&lat, r2)
  r2@units = "km"
  print(sprintf("%.1f", r2))
ip = 0
do d=0,ntmax-1
  if (ip .gt. plot_num-1) then
    break
  end if
  ind_t2 = ind(t.eq.tr(d+t_offset))
  if (d+t_offset.eq.0) then
    ind_t1 = 0
  else
    ind_t1 = ind(t.eq.tr(d+t_offset-1))
  end if
  u1   = tw1(ind_t2,:,:,:)
  u1   = dim_avg_n(tw1(ind_t1:ind_t2,:,:,:), 0)
  u2   = tw2(ind_t2,:,:,:)
  u2   = dim_avg_n(tw2(ind_t1:ind_t2,:,:,:), 0)
  u    = u1 
  ;u    = u2 - u1
  u    = u1 - u2
  hr1  = hrr1(d+t_offset,:,:,:)
  hr2  = hrr2(d+t_offset,:,:,:)
  hr   = hr1
  ;hr   = hr2 - hr1
  hr   = hr1 - hr2
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
  res1@cnLevelSpacingF      =  .5                  ; contour spacing
  res1@cnMinLevelValF       =  -4.5                 ; min level
  res1@cnMaxLevelValF       =  4.5                 ; max level
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
  res2@cnLevelSpacingF      =  2.                  ; contour spacing
  res2@cnMinLevelValF       =  -10.                 ; min level
  res2@cnMaxLevelValF       =  10.                 ; max level
  res2@gsnContourNegLineDashPattern = 2
  res2@gsnContourZeroLineThicknessF = 2.0
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
