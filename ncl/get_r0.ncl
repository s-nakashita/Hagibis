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
  f     = addfile ("../../netcdf/tigge/2019/jma/rotate/np_glb_2019100912_mean_ug.nc","r")
  wks   = gsn_open_wks ("png", "u850_lat_12Zmean" )          ; send graphics to PNG file
  yyyy = 2019
  mm   = 10
  dd   = 12
  hh   = 12
  u    = f->var(:,{850},:,:)                             ; get u for January
  printVarSummary(u)
  ndim = dimsizes(u)
  ntim = ndim(0)
  ;nlev = ndim(1)
  nlat = ndim(1)
  nlon = ndim(2)
  t    = f->time
  ;lev  = f->level
  lat  = f->lat
  ;units= t@units
  ;opt  = 0
  ;opt@calendar = t@calendar
  ;vt   = cd_inv_calendar(yyyy, mm, dd, hh, 0, 0, units, opt)
  ;print(vt)
  date  = cd_string(t, "%Y%N%D%H")
  print(date)
  id   = ind(date.eq."2019101012")
  print(id)
  plot_num = id/2 + 1
  plot = new(2,graphic)
  ;ft = tofloat(t)
  ;ft@units = units
  ;ft@calendar = t@calendar
;----------------------------------------------------------------
; get circular mean tangential wind
;----------------------------------------------------------------
  vtan = new((/ntim, nlat/),float)
  vtan!0 = "time"
  vtan!1 = "lat"
  vtan&time = t
  vtan&lat = lat
  vtan = dim_avg_n(u, 2)
  vtan@long_name = "Tangential mean wind"
;----------------------------------------------------------------
; calcurate radial derivative d(vtan)/dr
;----------------------------------------------------------------
  dvdr = vtan
  dvdr = center_finite_diff(vtan, dble2flt(lat), False, 0)
  dvdr = -dvdr
  dvdr@long_name = "Radial derivative of tangential mean wind"
  printMinMax(dvdr, 0)

;----------------------------------------------------------------
; draw plot
;----------------------------------------------------------------
data = new((/2, nlat/),float)
data!1 = "lat"
data&lat = lat
;do i=0, 0; plot_num-1
  d = id
  ;w = new((/nlev,nlat,nlon/),float)
  ;copy_VarCoords(u(0,:,:,:), w)
  ;w   = (/u(2*d,:,:,:)/)                             ; get u for January
  ;time  = t(2*d)
  time  = t(id)
  ;print(time)
  ;date  = cd_calendar(time, -3)
  ;date  = cd_string(time, "%Y%N%D%H")
  ;w&level@units = "hPa"
  ;w_zm = new((/nlev,nlat/),float)
  ;w_zm = w(:,:,0)
  ;printVarSummary(w_zm)
  ;w_zm = dim_avg_n(w, 2)
  ;w_zm@long_name = "Tangential wind"
  ;printVarSummary(w_zm)

  res                      = True                 ; plot mods desired

  res@gsnDraw	             = False
  res@gsnFrame	           = False
  ;res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  ;res@cnLevelSpacingF      =  1.                  ; contour spacing
  ;res@cnMinLevelValF       =  1.                 ; min level
  ;res@cnMaxLevelValF       =  8.                 ; max level
  ;res@cnMonoLineThickness  = False
  ;res@cnLineThicknesses    = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
  ;res@cnLineLabelsOn       = True                 ; turn on line labels
  ;res@cnFillOn             = True                 ; turn on color fill
  ;res@cnFillPalette        = "sunshine_9lev"             ; choose colormap
  ;res@tiYAxisString        = "pressure (hPa)"   
  ;res@gsnRightString       = "1.0e-4 /s"
  res@gsnLeftString       = ""
  res@gsnRightString       = "m/s"
  res@tiMainString         = date(id)
  
; Note: u is already on pressure levels. If this were model data, it
; would be necessary to interpolate from the hybrid coordinates to 
; pressure levels.
  ;vrt_zm = vrt_zm * 1.e04
  ;plotv  = gsn_csm_xy(wks, lat({70:}), vtan(2*d,{70:}), res )   ; place holder
  res@xyLineColors = (/"blue", "red"/)
  res@xyDashPatterns = (/0, 0/)
  res@gsnYRefLine  = (/8.0, 4.0e-6/)
  res@gsnYRefLineColor = (/"blue", "red"/)
  res@gsnYRefLineDashPattern = 1
  res@gsnYRefLineThicknessF = 3.0
  r = new(nlat, typeof(lat))
  deg2r(lat, r)
  print(r)
  rd = r(maxind(vtan(d, :)))
  print(rd)
  rd2 = 1.5*rd
  rf = r(minind(abs(vtan(d, :)-8.0)))
  r0 = 1.25*rf
  lat0 = r0
  r2deg(r0, lat0)
  print(lat0)
  res@gsnXRefLine = (/rd, rd2, rf, r0/)
  res@gsnXRefLineColor = (/"green", "green", "purple", "purple"/)
  res@gsnXRefLineDashPatterns = (/1, 0, 1, 0/)
  res@gsnXRefLineThicknessF = 3.0
  res@tiXAxisString = "radius(km)"
  res@tiYAxisString = "m/s"
  ;res@xyLabelMode = "Custom"
  ;res@xyExplicitLabels = (/"Vtan","dV/dr"/)
  ;plotdv  = gsn_csm_xy(wks, lat({70:}), dvdr(2*d,{70:}), res )   ; place holder
  ;overlay(plotv, plotdv)
  data(0, :) = vtan(d,:)
  data(1, :) = dvdr(d,:)
  plot(0) = gsn_csm_xy(wks, r(ind(lat.eq.70):), data(:,{70:}), res )
  E = new(nlat, typeof(r))
  calc_e(r, r0, E)
  rese                      = True                 ; plot mods desired

  rese@gsnDraw	            = False
  rese@gsnFrame	            = False
  rese@tiXAxisString        = "r/r0"
  rese@tiYAxisString        = "E(r)"
  rese@tiMainString         = date(id)
  imax = minind(abs(r-r0))
  plot(1) = gsn_csm_xy(wks, r(imax:)/r0, E(imax:), rese)
;end do  
  gsn_panel(wks, plot, (/1,2/), False)
  ;draw(plot(0))
  ;frame(wks)
end
