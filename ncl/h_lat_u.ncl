;----------------------------------------------------------------------
; h_lat_6.ncl
;----------------------------------------------------------------------
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
;----------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
  ;f     = addfile ("../../netcdf/gsm/gl/2019/10/rotate/np_init_vort.nc","r")
  f     = addfile ("../../netcdf/tigge/2019/jma/rotate/np_2019100900_mean.nc","r")
  wks   = gsn_open_wks ("png", "u_h" )          ; send graphics to PNG file
  yyyy = 2019
  mm   = 10
  dd   = 12
  hh   = 12
  ;u    = f->var(:,:,:,:)                             ; get u for January
  u    = f->UGRD(:,:,:,:)                             ; get u for January
  printVarSummary(u)
  ndim = dimsizes(u)
  ntim = ndim(0)
  nlev = ndim(1)
  nlat = ndim(2)
  nlon = ndim(3)
  t    = f->time
  lev  = f->level
  lat  = f->lat
  ;units= t@units
  ;opt  = 0
  ;opt@calendar = t@calendar
  ;vt   = cd_inv_calendar(yyyy, mm, dd, hh, 0, 0, units, opt)
  ;print(vt)
  date  = cd_string(t, "%Y%N%D%H")
  print(date)
  id   = ind(date.eq."2019101212")
  print(id)
  plot_num = id/2 + 1
  plot = new(plot_num,graphic)
  ;ft = tofloat(t)
  ;ft@units = units
  ;ft@calendar = t@calendar
  
do d=0,plot_num-1
  w = new((/nlev,nlat,nlon/),float)
  copy_VarCoords(u(0,:,:,:), w)
  w   = (/u(2*d,:,:,:)/)                             ; get u for January
  time  = t(2*d)
  print(time)
  ;date  = cd_calendar(time, -3)
  ;date  = cd_string(time, "%Y%N%D%H")
  ;w&level@units = "hPa"
  w_zm = new((/nlev,nlat/),float)
  w_zm = w(:,:,0)
  ;printVarSummary(w_zm)
  w_zm = dim_avg_n(w, 2)
  w_zm@long_name = "Tangential wind"
  printVarSummary(w_zm)

  res                      = True                 ; plot mods desired

  res@gsnDraw	             = False
  res@gsnFrame	           = False
  ;res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  ;res@cnLevelSpacingF      =  1.                  ; contour spacing
  ;res@cnMinLevelValF       =  1.                 ; min level
  ;res@cnMaxLevelValF       =  8.                 ; max level
  res@cnMonoLineThickness  = False
  res@cnLineThicknesses    = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
  res@cnLineLabelsOn       = True                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = "sunshine_9lev"             ; choose colormap
  res@tiYAxisString        = "pressure (hPa)"   
  ;res@gsnRightString       = "1.0e-4 /s"
  res@gsnLeftString       = ""
  res@gsnRightString       = "m/s"
  res@tiMainString         = date(2*d)
  
; Note: u is already on pressure levels. If this were model data, it
; would be necessary to interpolate from the hybrid coordinates to 
; pressure levels.
  ;vrt_zm = vrt_zm * 1.e04
  plot(d)  = gsn_csm_pres_hgt(wks, w_zm(:,{:80}), res )   ; place holder
end do  
  gsn_panel(wks, plot, (/2,4/), False)
end
