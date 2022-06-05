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
  init = "2019100900"
  fh1   = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/np_fcst_asia_"+init+".nc", "r")
  fh2   = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/np_fcst_asia_"+init+".nc","r")
  fh3   = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/np_fcst_asia_"+init+".nc", "r")
  yyyymmddhh = "2019101112"
  wks   = gsn_open_wks ("png", "vtan850_gsmtl959_i"+init+"-v"+yyyymmddhh )          ; send graphics to PNG file
  vtan1    = fh1->U(:,{850},:,:)                             ; get u for January
  vtan2    = fh2->U(:,{850},:,:)                             ; get u for January
  vtan3    = fh3->U(:,{850},:,:)                             ; get u for January
  ndim = dimsizes(vtan1)
  ntim = ndim(0)
  ;nlev = ndim(1)
  nlat = ndim(1)
  nlon = ndim(2)
  t    = fh1->time
  ;lev  = f->level
  lat  = fh1->lat
  ;units= t@units
  ;opt  = 0
  ;opt@calendar = t@calendar
  ;vt   = cd_inv_calendar(yyyy, mm, dd, hh, 0, 0, units, opt)
  ;print(vt)
  date  = cd_string(t, "%Y%N%D%H")
  print(date)
  id   = ind(date.eq.yyyymmddhh)
  print(id)
  plot = new(3, graphic)
  expname = (/"clim","est","mgd"/)
  var = (/vtan1,vtan2,vtan3/)
;----------------------------------------------------------------
; get circular mean tangential wind
;----------------------------------------------------------------
do i = 0, 2
  vtan = new((/ntim, nlat/),float)
  vtan!0 = "time"
  vtan!1 = "lat"
  vtan&time = t
  vtan&lat = lat
  vtan = dim_avg_n(var(i,:,:,:), 2)
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
  res@tiMainString         = expname(i)
  
; Note: u is already on pressure levels. If this were model data, it
; would be necessary to interpolate from the hybrid coordinates to 
; pressure levels.
  ;vrt_zm = vrt_zm * 1.e04
  ;plotv  = gsn_csm_xy(wks, lat({70:}), vtan(2*d,{70:}), res )   ; place holder
  res@xyLineColors = (/"blue", "red"/)
  res@xyDashPatterns = (/0, 0/)
  res@xyLineThicknessF = 3.0
  res@gsnYRefLine  = (/6.0, 4.0e-6/)
  res@gsnYRefLineColor = (/"blue", "red"/)
  res@gsnYRefLineDashPattern = 1
  res@gsnYRefLineThicknessF = 2.0
  r = new(nlat, typeof(lat))
  deg2r(lat, r)
  print(r)
  rd = r(maxind(vtan(d, :)))
  print(rd)
  rd2 = 1.5*rd
  rf = r(minind(abs(vtan(d, :)-6.0)))
  r0 = 1.25*rf
  lat0 = r0
  r2deg(r0, lat0)
  print(lat0)
  res@gsnXRefLine = (/rd, rd2, rf, r0/)
  res@gsnXRefLineColor = (/"green", "green", "purple", "purple"/)
  res@gsnXRefLineDashPatterns = (/1, 0, 1, 0/)
  res@gsnXRefLineThicknessF = 2.0
  res@tiXAxisString = "radius(km)"
  res@tiYAxisString = "m/s"
  ;res@xyLabelMode = "Custom"
  ;res@xyExplicitLabels = (/"Vtan","dV/dr"/)
  ;plotdv  = gsn_csm_xy(wks, lat({70:}), dvdr(2*d,{70:}), res )   ; place holder
  ;overlay(plotv, plotdv)
  data(0, :) = vtan(d,:)
  data(1, :) = dvdr(d,:)
  plot(i) = gsn_csm_xy(wks, r(:), data(:,:), res )
  ;E = new(nlat, typeof(r))
  ;calc_e(r, r0, E)
  ;rese                      = True                 ; plot mods desired

  ;rese@gsnDraw	            = False
  ;rese@gsnFrame	            = False
  ;rese@tiXAxisString        = "r/r0"
  ;rese@tiYAxisString        = "E(r)"
  ;rese@tiMainString         = date(id)
  ;imax = minind(abs(r-r0))
  ;plot(1) = gsn_csm_xy(wks, r(imax:)/r0, E(imax:), rese)
end do  
  gsn_panel(wks, plot, (/2,2/), False)
  ;draw(plot(0))
  ;frame(wks)
end
