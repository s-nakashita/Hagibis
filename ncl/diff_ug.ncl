;********************************************
; h_long_5.ncl
;********************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;********************************************
begin
  f00 = addfile("../../netcdf/tigge/2019/jma/rotate/inv_r8_glb_2019100900_mean_ug.nc", "r")
  f12 = addfile("../../netcdf/tigge/2019/jma/rotate/inv_r8_glb_2019100912_mean_ug.nc", "r")
  fa  = addfile("../../netcdf/gsm/gl/2019/10/rotate/inv_r8_init_ug.nc", "r")
  ;f     = addfile ("uvt.nc" , "r")
  lev = f00->level 
  leva= fa->level 
  lat = f00->lat
  lon = f00->lon
  t00 = f00->time
  t12 = f12->time
  ta  = fa->time
  date00 = cd_string(t00,"%Y%N%D%H")
  date12 = cd_string(t12,"%Y%N%D%H")
  datea  = cd_string(ta ,"%Y%N%D%H")
  d00    = ind(date00.eq."2019101012")
  d12    = ind(date12.eq."2019101012")
  da     = ind(datea.eq."2019101012")
  u00    = f00->var(d00,:,{20:25},{130:180})
  u12    = f12->var(d12,:,{20:25},{130:180})
  ua     = fa->var(da,:,{20:25},{130:180})
  
  u00_m = dim_avg_n(u00, 1)
  u00_m!0 = "lev"
  u00_m!1 = "lon"
  u00_m&lev = lev
  u00_m&lon = lon({130:180})
  u00_m@units = "m/s"
  u12_m = u00_m
  u12_m = dim_avg_n(u12, 1)

  ua_m = dim_avg_n(ua, 1)
  ua_m!0 = "lev"
  ua_m!1 = "lon"
  ua_m&lev = leva
  ua_m&lon = lon({130:180})
  ua_m@units = "m/s"
  
; ===========================
; color plot
; ===========================
  wks   = gsn_open_wks ("png", "u500_20N-25N_101012" )          ; open workstation
  ;plot = new(3,graphic)
  
  res                      = True                  ; plot mods desired
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@tiMainString         = "Zonal Wind 20N-25N"        ; title
  
  res@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
  res@cnLevelSpacingF      = 5.0                   ; contour interval
  res@cnMinLevelValF       = -50.                  ; min level
  res@cnMaxLevelValF       =  50.                  ; max level
  res@cnLineLabelsOn       = True                  ; turn on line labels
  res@cnFillOn             = True                  ; turn on color fill
  res@cnFillPalette        = "BlWhRe"              ; set color map
  
  res@xyLineThicknessF     = 3.0
  res@tiYAxisString       = "m/s"
  res@gsnYRefLine  = (/0.0/)
  res@gsnYRefLineColor = (/"black"/)
  res@gsnYRefLineDashPattern = 1
  res@gsnYRefLineThicknessF = 2.0
; note this data is already on pressure levels. If this were model data,
; it would be necessary to interpolate the data from hybrid coordinates
; to pressure levels.

  ;plot = gsn_csm_pres_hgt(wks,dv(:,{-10:20},{155}),res)
  ;res@gsnCenterString = "GSM analysis"
  res@xyLineColor = "Red"
  ;plot(0) = gsn_csm_pres_hgt(wks,ua_m({850:200},:),res)
  plot = gsn_csm_xy(wks,ua_m&lon, ua_m({500},:),res)
  ;res@gsnCenterString = "0000 UTC"
  res@xyLineColor = "Blue"
  ;plot(1) = gsn_csm_pres_hgt(wks,u00_m,res)
  plot00 = gsn_csm_xy(wks,u00_m&lon, u00_m({500},:),res)
  overlay(plot, plot00)
  ;res@gsnCenterString = "1200 UTC"
  res@xyLineColor = "Green"
  ;plot(2) = gsn_csm_pres_hgt(wks,u12_m,res)
  plot12 = gsn_csm_xy(wks,u12_m&lon, u12_m({500},:),res)
  overlay(plot, plot12)
  
  gres = True
  gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                            ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 5      ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  lineres = True
  lineres@lgLineColors = (/"red","blue","green"/) ; line colors
  lineres@lgLineThicknesses = 2.5                        ; line thicknesses
  lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

  textres = True
  textres@lgLabels = (/"GSM analysis","0000 UTC","1200 UTC"/)  ; legend labels (required)

  plot = simple_legend(wks,plot,gres,lineres,textres)
  
  draw(plot)
  frame(wks)
; the "{}" allows you to use coordinate subscripting to select a particular 
; latitude vice array indice of the latitude you want.
  ;resP = True
  ;gsn_panel(wks, plot, (/1,3/), resP)
end

