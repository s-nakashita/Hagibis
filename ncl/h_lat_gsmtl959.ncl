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
  Exp = "mgdsst"
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
  ug   = f->U(:,:,:,:)                             ; get u for January
  vg   = f->V(:,:,:,:)                             ; get u for January
  t    = f->time
  ntmax = dimsizes(t)
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
  plot_num = 15
  plot = new(plot_num,graphic)
  ;ft = tofloat(t)
  ;ft@units = units
  ;ft@calendar = t@calendar
  date  = cd_string(t, "%Y%N%D%H")
  print(date)
  t_offset = 15
  wks   = gsn_open_wks ("png", "vor_h-lat_"+Exp+"_"+date(t_offset)+"-"+date(min((/t_offset+14,ntmax-1/))) )          ; send graphics to PNG file
  vor = ug
  vor = uv2vr_cfd(ug, vg, ug&lat, ug&lon, 3)
  vor@standard_name = "Vorticity"
  vor@units = "/s"
  r = new(dimsizes(ug&lat), typeof(ug&lat))
  deg2r(ug&lat, r)
  r@units = "km"
  print(r)
  print(sprintf("%.1f", r))
do d=0,plot_num-1
  if (d+t_offset .gt. ntmax-1) then
    break
  end if
  u   = vor(d+t_offset,:,:,:)                             ; get u for January
  time  = t(d+t_offset)
  print(time)
  ;date  = cd_calendar(time, -3)
  ;date  = cd_string(time, "%Y%N%D%H")
  u&level@units = "hPa"
  u_zm = u(:,:,0)
  printVarSummary(u_zm)
  u_zm = zonalAve(u)
  printVarSummary(u_zm)

  res                      = True                 ; plot mods desired

  res@gsnDraw	             = False
  res@gsnFrame	           = False
  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res@cnLevelSpacingF      =  1.                  ; contour spacing
  res@cnMinLevelValF       =  1.                 ; min level
  res@cnMaxLevelValF       =  12.                 ; max level
  res@cnMonoLineThickness  = False
  res@cnLineThicknesses    = (/1.0,1.0,1.0,1.0,3.0,1.0,1.0,1.0,1.0,3.0,1.0,1.0,1.0,1.0,3.0/)
  res@cnMonoLineDashPattern= False
  res@cnLineDashPatterns   = (/  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
  res@cnLineLabelsOn       = True                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = "WhiteYellowOrangeRed"             ; choose colormap
  res@tiYAxisString        = "pressure (" + u&level@units + ")"   
  res@tmYROn               = False
  res@tmYRLabelsOn         = False
  res@tmYUseLeft           = False
  res@trYReverse           = True
  ;res@tmXBMode             = "Explicit"
  ;res@tmXBValues           = r
  ;res@tmXBLabels           = r
  res@gsnRightString       = "1.0e-4 /s";u@units
  res@tiMainString         = date(d+t_offset)
  
; Note: u is already on pressure levels. If this were model data, it
; would be necessary to interpolate from the hybrid coordinates to 
; pressure levels.
  data = u_zm * 1.0e4
  data!0 = "level"
  data&level = u&level
  data!1 = "distance"
  data&distance = r
  data@standard_name = u_zm@standard_name
  ;plot(d)  = gsn_csm_pres_hgt(wks, u_zm({100:1000},{75:90}), res )   ; place holder
  plot(d)  = gsn_csm_contour(wks, data({100:1000},{0:500}), res )   ; place holder
end do  
  gsn_panel(wks, plot, (/3,5/), False)
end
