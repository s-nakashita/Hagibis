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
  f     = addfile ("../../netcdf/tigge/2019/jma/rotate/np_100900_mean_v.nc","r")
  wks   = gsn_open_wks ("png", "vor_h_lat_00Zmean" )          ; send graphics to PNG file
  yyyy = 2019
  mm   = 10
  dd   = 12
  hh   = 12
  v    = f->vort(:,:,:,:)                             ; get u for January
  t    = f->time
  units= t@units
  opt  = 0
  opt@calendar = t@calendar
  vt   = cd_inv_calendar(yyyy, mm, dd, hh, 0, 0, units, opt)
  print(vt)
  id   = ind(t.eq.vt)
  print(id)
  plot_num = id/2 + 1
  plot = new(plot_num,graphic)
  ft = tofloat(t)
  ft@units = units
  ft@calendar = t@calendar
  date  = cd_string(ft, "%Y%N%D%H")
  print(date)
do d=0,plot_num-1
  vrt   = v(2*d,:,:,:)                             ; get u for January
  time  = t(2*d)
  print(time)
  ;date  = cd_calendar(time, -3)
  ;date  = cd_string(time, "%Y%N%D%H")
  vrt&level@units = "hPa"
  vrt_zm = vrt(:,:,0)
  printVarSummary(vrt_zm)
  vrt_zm = zonalAve(vrt)
  printVarSummary(vrt_zm)

  res                      = True                 ; plot mods desired

  res@gsnDraw	             = False
  res@gsnFrame	           = False
  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res@cnLevelSpacingF      =  1.                  ; contour spacing
  res@cnMinLevelValF       =  1.                 ; min level
  res@cnMaxLevelValF       =  8.                 ; max level
  res@cnMonoLineThickness  = False
  res@cnLineThicknesses    = (/1.0,1.0,1.0,1.0,3.0,1.0,1.0,1.0/)
  res@cnLineLabelsOn       = True                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = "sunshine_9lev"             ; choose colormap
  res@tiYAxisString        = "pressure (" + vrt&level@units + ")"   
  res@gsnRightString       = "1.0e-4 /s"
  res@tiMainString         = date(2*d)
  
; Note: u is already on pressure levels. If this were model data, it
; would be necessary to interpolate from the hybrid coordinates to 
; pressure levels.
  vrt_zm = vrt_zm * 1.e04
  plot(d)  = gsn_csm_pres_hgt(wks, vrt_zm(:,:20), res )   ; place holder
end do  
  gsn_panel(wks, plot, (/2,4/), False)
end
