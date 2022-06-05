;----------------------------------------------------------------------
; h_lat_7.ncl
;----------------------------------------------------------------------
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
;----------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
  level = "surf"
;---Read in variables
  if (level .eq. "surf") then
    f    = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/fcst_surf_asia_2019100912.nc","r")
    fp   = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p/fcst_surf_asia_2019100912.nc","r")
    ps   = f->PSEA(:,{10:15},:)
    psp  = fp->PSEA(:,{10:15},:)
    zp = psp
    zp = psp - ps
  else
    f    = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/fcst_p_asia_2019100912.nc","r")
    fp   = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p/fcst_p_asia_2019100912.nc","r")
    z    = f->Z(:,{level},{10:15},:)
    zp   = fp->Z(:,{level},{10:15},:)
    zp = zp - z ;diff
  end if
  time = f->time
  date = cd_string(time, "%N%D %H:%M")
  tmax = cd_inv_calendar(2019,10,11,6,0,0,time@units,0)
  indmax = ind(time.eq.tmax)
  print(indmax)
;*************************************
; averaged between 10N and 15N
;*************************************
  zd = dim_avg_n(zp,1)
  zd!0 = "time"
  zd!1 = "lon"
  zd&time = time
  zd&lon = zp&lon
  ;zd&time@units = time@units
  zd&lon@units = "degrees_east"
  if (level.eq."surf") then
    zd@standard_name = "Sea level pressure"
  else
    zd@standard_name = "Geopotential height"
  end if
;***************************************
; PLOT
;**************************************
  if (level.eq."surf") then
    pltName = "t_lon_diff_surf"
  else
    pltName = "t_lon_diff_"+sprinti("%3d", level)+"hPa"
  end if
  wks   = gsn_open_wks ("png", pltName )          ; send graphics to PNG file
  
  cmap = read_colormap_file("GMT_polar")
  ;nc   = dimsizes(cmap(:,0))
  res                      = True
  ;res@gsnDraw              = False
  ;res@gsnFrame             = False
  ;res@tiMainString         = zd@standard_name
  res@cnFillPalette        = cmap;(10:nc-3,:)
  res@cnFillOn             = True                ; color on
  res@cnLevelSelectionMode = "ManualLevels"      ; manual contouring
  if (level .eq. "surf") then
  res@cnMinLevelValF       = -2.0               ; set min contour level
  res@cnMaxLevelValF       =  2.0               ; set max contour level
  res@cnLevelSpacingF      = .2                 ; set contour spacing
  else
  res@cnMinLevelValF       = -10.0               ; set min contour level
  res@cnMaxLevelValF       =  10.0               ; set max contour level
  res@cnLevelSpacingF      = 1.                  ; set contour spacing
  end if
  res@lbLabelStride        = 2                   ; every other label bar label
  res@tmYLOn = True
  res@tmYLMode = "Explicit"
  ;res@tmYLLabelFontHeightF = 0.005
  if (level .eq. "surf") then
  res@tmYLValues = time(:indmax:6)
  res@tmYLLabels = date(:indmax:6)
  else
  res@tmYLValues = time(:indmax:2)
  res@tmYLLabels = date(:indmax:2)
  end if
  plot  = gsn_csm_hov(wks, zd(:indmax,{120:180}), res )
  ;res2                     = True
  ;;res2@gsnDraw             = False
  ;res2@gsnFrame            = False
  ;res2@cnFillOn            = False
  ;res2@cnLinesOn           = True
  ;res2@cnLineThicknessF    = 3.0
  ;res2@cnLevelSelectionMode= "ManualLevels"      ; manual contouring
  ;res2@cnMinLevelValF      = -30.0               ; set min contour level
  ;res2@cnMaxLevelValF      =  30.0               ; set max contour level
  ;res2@cnLevelSpacingF     = 5.                  ; set contour spacing
  ;res2@gsnContourNegLineDashPattern = 2
  ;res2@gsnContourZeroLineThicknessF = 5.0
  ;res2@gsnRightString = up@standard_name
  ;res2@gsnLeftString = ""

  ;plot  = gsn_csm_pres_hgt(wks, zd({:100},{120:180}), res )
  ;plot2 = gsn_csm_pres_hgt(wks, up({:100},{120:180}), res2 )
  ;overlay(plot, plot2)
  ;frame(wks)
end
