load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "calc_e.ncl"
begin
 ;rdrdir = "/Volumes/dandelion/netcdf/JMAGPV/hourly/2019/10"
 gpvdir = "/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est"
 init = "2019100912"

 ndate = 25
 ;ndate = 29
 time = new(ndate, float)
 time@units = "seconds since 2019-10-09 12:00:00"
 do it = 0, ndate-1
   time(it) = it*3.*60.*60.
 end do
; GSM GPV
 files = (/gpvdir+"/np_fcst_asia_"+init+".nc", \ ; CNTL
           gpvdir+"_0912+p/np_fcst_asia_"+init+".nc", \ ; EN-
           gpvdir+"_0912+p2/np_fcst_asia_"+init+".nc",\ ; TC
           gpvdir+"_0912+pn/np_fcst_asia_"+init+".nc"/) ; EN+
 ltxts = (/"wmax_"+init+"_cntl.txt", \
           "wmax_"+init+"_en-.txt", \
           "wmax_"+init+"_tc.txt", \
           "wmax_"+init+"_en+.txt"/)
 rtxts = (/"rmw_"+init+"_cntl.txt", \
           "rmw_"+init+"_en-.txt", \
           "rmw_"+init+"_tc.txt", \
           "rmw_"+init+"_en+.txt"/)
 r25txts = (/"r25_"+init+"_cntl.txt", \
             "r25_"+init+"_en-.txt", \
             "r25_"+init+"_tc.txt", \
             "r25_"+init+"_en+.txt"/)
 v_thres = 25.0 ; strong wind threshold
 do j = 0, dimsizes(files)-1
  f = addfile(files(j), "r")
  ;fr = addfile(rfiles(j), "r")
  t = f->time
  lat = dble2flt(f->lat)
  lon = f->lon
  nlat = dimsizes(lat)
  nlon=dimsizes(lon)
  dnlon = nlon / 8
  r = new(nlat, typeof(lat))
  deg2r(lat, r)
  print(sprintf("%.1f",r))
  r@units = "km"
  ldata = new(ndate, float)
  rdata = new(ndate, typeof(r))
  r25data = new(ndate, typeof(r))
  do it = 0, ndate-1
    vmax = 0.0
    rmax = 0.0
    imax = 0
    var = f->U(it,{850},:,:)
    vazim = var(:,0)
    vazim = dim_avg_n(var,1)
    vmax = max(vazim)
    imax = maxind(vazim)
    rmax = r(imax)
    do ilat = 1, nlat-1
      if (vazim(ilat) .lt. v_thres) then
        break
      end if
    end do
    wgt = (v_thres - vazim(ilat))/(vazim(ilat-1)-vazim(ilat))
    r25 = r(ilat-1)*wgt + r(ilat)*(1.0-wgt)
;    rtmp = new(8, typeof(r))
;    do i = 0, 7
;      ;var = f->U10m(it,:,i*dnlon:(i+1)*dnlon-1)
;      var = f->U(it,{850},:,i*dnlon:(i+1)*dnlon-1)
;      tmp = dim_avg_n(var,1)
;      vtmp = max(tmp)
;      itmp = maxind(tmp)
;      rtmp = r(itmp)
;      if (vtmp .gt. vmax) then
;        vmax = vtmp
;      ;  rmax = rtmp
;        imax = itmp
;      end if
;      do ilat = 1, nlat-1
;        if (tmp(ilat) .lt. v_thres) then
;          break
;        end if
;      end do
;      wgt = (v_thres - tmp(ilat))/(tmp(ilat-1)-tmp(ilat))
;      rtmp(i) = r(ilat-1)*wgt + r(ilat)*(1.0-wgt)
;    end do
;    rmax = max(rtmp)
    ldata(it) = vmax
    rdata(it) = rmax
    r25data(it) = r25
    print(imax)
    print(rmax)
    print(r25)
  end do
  ldata!0 = "time"
  ldata&time = time
  printMinMax(ldata, 0)
  asciiwrite(ltxts(j), ldata)
  rdata!0 = "time"
  rdata&time = time
  printMinMax(rdata, 0)
  asciiwrite(rtxts(j), rdata)
  r25data!0 = "time"
  r25data&time = time
  printMinMax(r25data, 0)
  asciiwrite(r25txts(j), r25data)
  delete(ldata)
  delete(rdata)
  delete(r25data)
  delete(t)
 end do
end