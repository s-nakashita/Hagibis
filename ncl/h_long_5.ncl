;********************************************
; h_long_5.ncl
;********************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;********************************************
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
  wks   = gsn_open_wks ("png", "h_lon_merra2_20191006-12" )          ; open workstation
  plot = new(8,graphic)
  res                      = True                  ; plot mods desired
  res@gsnFrame = False
  res@gsnDraw  = False

  res@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
  res@cnLevelSpacingF      = 2.5                   ; contour interval
  res@cnMinLevelValF       = -25.                  ; min level
  res@cnMaxLevelValF       =  25.                  ; max level
  res@cnLineLabelsOn       = True                  ; turn on line labels
  res@cnFillOn             = True                  ; turn on color fill
  res@cnFillPalette        = "BlWhRe"              ; set color map
  d = 0
  do dd=6,12
  yyyymmdd = "201910"+sprinti("%0.2d",dd)
  f = addfile("/Volumes/dandelion/netcdf/nasa/merra2/MERRA2_400.tavg3_3d_asm_Nv."+yyyymmdd+".nc4","r")
  ;f     = addfile ("uvt.nc" , "r") 
  time = f->time
  date = cd_string(time, "%D")
  lat = f->lat
  lon = f->lon 
  u1     = f->U(:,:,:,:)
  ;u2     = f->U(:,:,:,:)
  pl1    = f->PL(:,:,:,:)
  ;pl2    = f->PL(:,:,:,:)
  ;u    = u1
  ;u    = 0.5*(u1+u2)
  ;pl   = pl1
  ;pl   = 0.5*(pl1+pl2)
  ;print(u1&lev)
  pnew = (/10000.0,15000.0,20000.0,25000.0,30000.0,40000.0,50000.0,60000.0,70000.0,85000.0,90000.0,92500.0,95000.0,100000.0/)
  u  = int2p_n(pl1, u1, pnew, 0, 1)
  ;up2  = int2p_n(pl2, u2, pnew, 0, 1)
  ;u    = up1
  ;u    = 0.5*(up1+up2)
  pnew = pnew*0.01 ;hPa
  pnew@units = "hPa"
  u!0 = "time"
  u!1 = "lev"
  u!2 = "lat"
  u!3 = "lon"
  u&time = time
  u&lev = pnew
  u&lat = lat
  u&lon = lon
  ;vp1  = int2p_n(pl1, v1, pnew, 1, 0)
  ;vp2  = int2p_n(pl2, v2, pnew, 1, 0)
  ;v    = vp1
  ;v    = 0.5*(vp1+vp2)

  ;dv = u
  ;dv = uv2dvF(u, v) ;divergence
  ;dv = uv2vrF(u, v) ;vorticity
  ;sfvp = uv2sfvpF(u, v)
  ;dv = (/sfvp(1,:,:,:)/) ; velocity potential
  ;dvp  = int2p_n(pl(:,{-10:20},{150:170}), dv(:,{-10:20},{150:170}), pnew, 1, 0)
  ;dvp = dvp * 1e6
  ;dvp!2 = "lon"
  ;dvp&lon = lon
  ;dv_z = dim_avg_n(dvp, 2)
  ;dv_z!0 = "lev"
  ;dv_z!1 = "lat"
  ;;dv!2 = "lon"
  ;dv_z&lev = pnew
  ;dv_z&lat = lat({-10:20})
  ;;dv&lon = lon
  ;dv_z&lev@units = "Pa"
  ;dv_z&lat@units = "degrees_north"
  um = dim_avg_n(u(:,:,{-10:20},:),2)
  u_z = dim_avg_n(um,0)
  u_z!0 = "level"
  u_z!1 = "lon"
  u_z&level = pnew
  u_z&lon = lon 
  if (d.eq.0) then
   u_zm = u_z
  else 
   u_zm = u_zm + u_z
  end if
; ===========================
; color plot
; ===========================
  
; note this data is already on pressure levels. If this were model data,
; it would be necessary to interpolate the data from hybrid coordinates
; to pressure levels.

  ;plot = gsn_csm_pres_hgt(wks,dv(:,{-10:20},{155}),res)
  res@gsnCenterString = date(0)
  plot(d) = gsn_csm_pres_hgt(wks,u_z(:,{110:170}),res)
  d = d + 1
  end do
  res@gsnCenterString = "time average"
  u_zm = u_zm / d 
  plot(d) = gsn_csm_pres_hgt(wks,u_zm(:,{110:170}),res)
  resP = True
  resP@gsnPanelMainString         = "U average 10S-20N"        ; title
  gsn_panel(wks,plot,(/2,4/),resP)
end

