;*************************************************
; polar_8.ncl
;************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;************************************************

;************************************************
; Function to retrieve the location of a plot so we can calcurate
; a location for a title on the left side of the plot.
;************************************************
function get_title_position(plot)
local vpy, vpx, vph
begin
  getvalues plot
      "vpYF"		: vpy
      "vpXF"		: vpx
      "vpHeightF"	: vph
  end getvalues
  pos = (/vpx*0.75,vpy-(vph/2.)/)
  return(pos)
end

;************************************************
; Main code
;************************************************
begin
;************************************************
; read in atmospheric data
;************************************************
  yyyymmddhh = "2019100912"
  yyyy 	     = str_get_cols(yyyymmddhh,0,3)
  mm  	     = str_get_cols(yyyymmddhh,4,5)
  dd  	     = str_get_cols(yyyymmddhh,6,7)
  hh  	     = str_get_cols(yyyymmddhh,8,9)
  ;Exp = "clim"
  if (Exp .eq. "cntl") then
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/np_fcst_asia_"+yyyymmddhh+".nc"
  else if (Exp .eq. "p") then
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p/np_fcst_asia_"+yyyymmddhh+".nc"
  else if (Exp .eq. "p2") then
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p2/np_fcst_asia_"+yyyymmddhh+".nc"
  else if (Exp .eq. "pn") then
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+pn/np_fcst_asia_"+yyyymmddhh+".nc"
  else if (Exp .eq. "pf") then
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+pf/np_fcst_asia_"+yyyymmddhh+".nc"
  end if
  end if
  end if
  end if
  end if
  a = addfile(f,"r")
  t = a->time
  ntime = dimsizes(t)
  lon = a->lon
  nlon = dimsizes(lon)
  print(nlon)
  ;;       850    800    700    600    500    400    300
  ;wgt = (/1./22.,3./22.,2./11.,2./11.,2./11.,2./11.,1./11./)
  ;print(wgt)
  lev = a->level({1000:100});({850:300})
  print(lev)
  nlev = dimsizes(lev)
  print(nlev)
  ;dp   = dpres_plevel_Wrap(lev, lev({850}), lev({300}), 0)
  ;dpall= lev({850}) - lev({300})
  ;wgt  = dble2flt(dp / dpall)
  ;print(wgt)

  odir = "steer_" + Exp + "_" + yyyymmddhh
  system("mkdir -p " + odir)
  i = 0
  do d = 0, ntime-1
    date  = cd_string(t(d), "%Y%N%D%H")
    print(date)
    if (Exp.eq."era5") then
    ;u = a->u(d,{850:300},:,:)
    ;v = a->v(d,{850:300},:,:)
    u = a->u(d,{1000:100},:,:)
    v = a->v(d,{1000:100},:,:)
    else
    ;u = a->U(d,{850:300},:,:)
    ;v = a->V(d,{850:300},:,:)
    u = a->U(d,{1000:100},:,:)
    v = a->V(d,{1000:100},:,:)
    end if
    assignFillValue(v,u)                           ; in contributed.ncl   
    assignFillValue(u,v)
    uzav = dim_avg_n(u,2) ;axi-symmetrical value
    vzav = dim_avg_n(v,2)
    ;subtract axi-symmetrical components
    u_ = u 
    u_ = u_ - conform_dims(dimsizes(u), uzav, (/0,1/)) ; level x lat x lon
    v_ = v
    v_ = v_ - conform_dims(dimsizes(v), vzav, (/0,1/)) ; level x lat x lon
    ;;vertical weighted summation
    ;u_ = u_ * conform_dims(dimsizes(u_), wgt, 0)
    ;v_ = v_ * conform_dims(dimsizes(v_), wgt, 0)
    ;upav = dim_sum_n(u_, 0)
    ;upav!0 = "lat"
    ;upav!1 = "lon"
    ;upav&lat = u&lat
    ;upav&lon = u&lon
    ;vpav = upav
    ;vpav = dim_sum_n(v_, 0)
    ;averaged between 3 to 5 degree of latitudes
    ;us = dim_avg_n(upav({87:85},:), 0)
    ;us!0 = "lon"
    ;us&lon = u&lon
    us = dim_avg_n(u_(:,{87:85},:), 1)
    us!0 = "level"
    us&level = u&level
    us!1 = "lon"
    us&lon = u&lon
    us@long_name = "Tangential component of steering flow"
    us@units     = "m/s"
    ;vs = dim_avg_n(vpav({87:85},:), 0)
    ;vs!0 = "lon"
    ;vs&lon = u&lon
    vs = dim_avg_n(v_(:,{87:85},:), 1)
    vs!0 = "level"
    vs&level = v&level
    vs!1 = "lon"
    vs&lon = v&lon
    vs@long_name = "Radial component of steering flow"
    vs@units     = "m/s"
    ;asciiwrite(odir+"/us_"+date+".txt", us)
    ;asciiwrite(odir+"/vs_"+date+".txt", vs)
    filename = odir+"/ushod_"+date+".mtx"
    fmt_str = "%s" + str_concat(conform_dims(nlev,"%8.3f",-1))
    row_labels = sprinti("%4i",ispan(1,nlon,1))
    dlist = NewList("lifo")
    do j=nlev-1,0,1
      ListPush(dlist,(/lev(j)/))
    end do
    ListPush(dlist,"level")
    write_table(filename, "w", dlist, fmt_str)
    delete(dlist)
    dlist = NewList("lifo")
    do j=nlev-1,0,1
      ListPush(dlist,(/us(j,:)/))
    end do
    ListPush(dlist,row_labels)
    write_table(filename, "a", dlist, fmt_str)
    delete(dlist)
    filename = odir+"/vshod_"+date+".mtx"
    fmt_str = "%s" + str_concat(conform_dims(nlev,"%8.3f",-1))
    row_labels = sprinti("%4i",ispan(1,nlon,1))
    dlist = NewList("lifo")
    do j=nlev-1,0,1
      ListPush(dlist,(/lev(j)/))
    end do
    ListPush(dlist,"level")
    write_table(filename, "w", dlist, fmt_str)
    delete(dlist)
    dlist = NewList("lifo")
    do j=nlev-1,0,1
      ListPush(dlist,(/vs(j,:)/))
    end do
    ListPush(dlist,row_labels)
    write_table(filename, "a", dlist, fmt_str)
  end do
end
