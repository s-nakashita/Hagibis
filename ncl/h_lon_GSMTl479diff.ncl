;----------------------------------------------------------------------
; h_lat_7.ncl
;----------------------------------------------------------------------
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
;----------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin

;---Read in variables
  f    = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/fcst_p_asia_2019100912.nc","r")
  fp   = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p/fcst_p_asia_2019100912.nc","r")
  ;d = 0 ; 3h after
  z    = f->Z(d,:,{10:15},:)
  zp   = fp->Z(d,:,{10:15},:)
  u    = fp->U(d,:,{10:15},:)
  lev  = f->level          ; pressure levels (millibars)
  time = f->time
  date = cd_string(time, "%D%H")
  zp = zp - z ;diff
;*************************************
; averaged between 10N and 15N
;*************************************
  zd = dim_avg_n(zp,1)
  zd!0 = "level"
  zd!1 = "lon"
  zd&level = lev
  zd&lon = zp&lon
  zd&level@units = "hPa"
  zd&lon@units = "degrees_east"
  zd@standard_name = "Geopotential height"
  up = zd
  up = dim_avg_n(u,1)
  up@standard_name = "Zonal wind"
  vp = up
  vp = 0.
;***************************************
; PLOT
;**************************************
  wks   = gsn_open_wks ("png", "h_lon_zdiff+u_"+date(d) )          ; send graphics to PNG file
  
  cmap = read_colormap_file("GMT_polar")
  ;nc   = dimsizes(cmap(:,0))
  res                      = True
  ;res@gsnDraw              = False
  ;res@gsnFrame             = False
  res@tiMainString         = cd_string(time(d), "%N%D %H:%M")
  res@cnFillPalette        = cmap;(10:nc-3,:)
  res@cnFillOn             = True                ; color on
  res@cnLevelSelectionMode = "ManualLevels"      ; manual contouring
  res@cnMinLevelValF       = -20.0               ; set min contour level
  res@cnMaxLevelValF       =  20.0               ; set max contour level
  res@cnLevelSpacingF      = 2.                  ; set contour spacing
  res@lbLabelStride        = 2                   ; every other label bar label
  
  res@vcRefAnnoOrthogonalPosF  = -1.165          ; move reference vec up
  res@vcRefAnnoParallelPosF    = 0.995
  res@pmLabelBarOrthogonalPosF = -0.05           ; move color label bar up
  res@vcRefMagnitudeF          = 30.0           ; define vector ref mag
  res@vcRefLengthF             = 0.045          ; define length of vec ref
  res@vcMinDistanceF           = 0.03
  res@vcMinMagnitudeF          = 5.0
  res@vcLineArrowThicknessF    = 3.0
  res@vcMapDirection           = False

  plot  = gsn_csm_pres_hgt_vector(wks, zd({:100},{120:180}) \
                                     , up({:100},{120:180}) \
                                     , vp({:100},{120:180}) \
                                     , res )
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
