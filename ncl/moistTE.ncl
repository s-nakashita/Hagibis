;**********************************************************************
; Simple function that takes the difference of two variables and
; copies the coordinate arrays from one of the variables to the
; difference variable for plotting later.
;**********************************************************************
function diff_vars(x1,x2)
local xdiff
begin
  xdiff = x1-x2
  copy_VarCoords(x1,xdiff)
  return(xdiff)
end

;**********************************************************************
; Main code
;**********************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin 
  ;Reading Historical data

  fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/fcst_p_asia_2019100912.nc", "r")
  fh2     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p/fcst_p_asia_2019100912.nc","r")
  ;fs1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/fcst_surf_asia_2019100912.nc", "r")
  ;fs2     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p/fcst_surf_asia_2019100912.nc","r")
  ;fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_surf_2019100900.nc", "r")
   
  ;constant
  pi = acos(-1.0)
  deg2rad = pi / 180.0
  cp = 1005.7
  rd = 287.04
  lh = 2.5104e6
  tr = 270.0
  pr = 1000.0 ;hPa

  ;verification region
  lonmin = 137.0
  lonmax = 142.0
  latmin = 33.0
  latmax = 37.0
  ;dimensions
  level = fh1->level({1000:200})
  nlev = dimsizes(level)
  lat = fh1->lat({latmin:latmax})
  nlat = dimsizes(lat)
  lon = fh1->lon({lonmin:lonmax})
  nlon = dimsizes(lon)
  t3h = fh1->time
  ntim = dimsizes(t3h)
  hour = cd_convert(t3h, "hours since 2019-10-09 12:00")
  ;variables
  u   = fh1->U(:, {1000:200}, {latmin:latmax}, {lonmin:lonmax})
  up  = fh2->U(:, {1000:200}, {latmin:latmax}, {lonmin:lonmax})
  v   = fh1->V(:, {1000:200}, {latmin:latmax}, {lonmin:lonmax})
  vp  = fh2->V(:, {1000:200}, {latmin:latmax}, {lonmin:lonmax})
  t   = fh1->T(:, {1000:200}, {latmin:latmax}, {lonmin:lonmax})
  tp  = fh2->T(:, {1000:200}, {latmin:latmax}, {lonmin:lonmax})
  rh  = fh1->RH(:, {1000:200}, {latmin:latmax}, {lonmin:lonmax})
  rhp = fh2->RH(:, {1000:200}, {latmin:latmax}, {lonmin:lonmax})
  ps  = fh1->P(:, {latmin:latmax}, {lonmin:lonmax})
  psp = fh2->P(:, {latmin:latmax}, {lonmin:lonmax})
  print(level)
  print(lat)
  print(lon)
;****************************************************************************
; Calculate specific humidity (kg/kg)
;****************************************************************************
  q   = rh(:,:,:,:)
  q   = mixhum_ptrh(conform(t,level,1), t, rh, 2)
  qp  = rhp(:,:,:,:)
  qp  = mixhum_ptrh(conform(tp,level,1), tp, rhp, 2)
;****************************************************************************
; Difference
;****************************************************************************
  ud  = diff_vars(up,u)
  vd  = diff_vars(vp,v)
  td  = diff_vars(tp,t)
  qd  = diff_vars(qp,q)
  psd = diff_vars(psp,ps)
  printMinMax(ud,0)
  printMinMax(vd,0)
  printMinMax(td,0)
  printMinMax(qd,0)
  printMinMax(psd,0)
;****************************************************************************
; Calculate layer thickness factor
;****************************************************************************
  wlev = new(nlev, float)
  wlev = 0.0
  do ilev = 0, nlev-2
    levc = (level(ilev)+level(ilev+1))*0.5
    wlev(ilev) = wlev(ilev) + (level(ilev)-levc)/800.0
    wlev(ilev+1) = wlev(ilev+1) + (levc-level(ilev+1))/800.0
  end do 
  print(wlev)
;****************************************************************************
; Calculate area factor
;****************************************************************************
  wlat = new(nlat, float)
  dlat = lat(1) - lat(0)
  do ilat = 0, nlat-1
    wlat(ilat) = sqrt(cos(lat(ilat)*deg2rad))
  end do
  print(wlat)
;****************************************************************************
; Calculate other factors
;****************************************************************************
  wt  = sqrt(cp/tr)
  wps = sqrt(rd*tr)/pr
  wq  = lh/sqrt(cp*tr)
;****************************************************************************
; Multiply factor
;****************************************************************************
  ud  = ud  * conform(ud,wlev,1)
  ud  = ud  * conform(ud,wlat,2)
  vd  = vd  * conform(vd,wlev,1)
  vd  = vd  * conform(vd,wlat,2)
  td  = td  * conform(td,wlev,1)
  td  = td  * conform(td,wlat,2)
  td  = td  * conform(td,wt,-1)
  qd  = qd  * conform(qd,wlev,1)
  qd  = qd  * conform(qd,wlat,2)
  qd  = qd  * conform(qd,wq,-1)
  psd = psd * conform(psd,wlat,1)
  psd = psd * conform(psd,wps,-1)
  printMinMax(ud,0)
  printMinMax(vd,0)
  printMinMax(td,0)
  printMinMax(qd,0)
  printMinMax(psd,0)
;****************************************************************************
; Calculate moist total energy
;****************************************************************************
  te = new(ntim, float)
  te = 0.0
  te = te + dim_sum_n(ud*ud, (/1,2,3/))
  te = te + dim_sum_n(vd*vd, (/1,2,3/))
  te = te + dim_sum_n(td*td, (/1,2,3/))
  te = te + dim_sum_n(qd*qd, (/1,2,3/))
  te = te + dim_sum_n(psd*psd, (/1,2/))
  te = te*0.5
  write_table("TE_p.txt","w",[/hour,te/], \
  "%5.1f %16.8e")
end
