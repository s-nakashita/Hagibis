load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "calc_e.ncl"
begin
 ;rdrdir = "/Volumes/dandelion/netcdf/JMAGPV/hourly/2019/10"
 gpvdir = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"
 init = "2019100900"

 wks = gsn_open_wks("png", "twspd_avr_core_"+init)
 rmax = 200.0

 res = True
 res@gsnDraw = False
 res@gsnFrame = False
 res@gsnCenterString = ""
 res@gsnRightString  = ""
 res@gsnLeftString   = ""
 res@xyLineThicknessF = 4.0

 ndate = 15
 ;ndate = 29
 time = new(ndate, float)
 time@units = "seconds since 2019-10-09 00:00:00"
 do t = 0, ndate-1
   time(t) = t*6.*60.*60.
   ;time(t) = t*3.*60.*60.
 end do
 res@tmXBMode = "Explicit"
 ;res@tmXBValues = time(::4)
 ;res@tmXBLabels = cd_string(time(::4), "%D-%HZ")
 res@tmXBValues = time(::8)
 res@tmXBLabels = cd_string(time(::8), "%D-%HZ")
 res@tmXBMinorOn = True
 res@tmXBMinorValues = time
 res@trXMinF = time(0)
 res@trXMaxF = time(ndate-1)
  
 ; left:latent heat flux
 resL = res
 resL@trYMinF = 0.0
 resL@trYMaxF = 800.0
 resL@tiYAxisString = "Latent heat flux (W/m~S~2~N~)"
 
 ; right:wind speed
 resR = res
 resR@xyDashPatterns = 1
 resR@trYMinF = 15.0
 resR@trYMaxF = 35.0
 ;res@gsnRightString = ""
 ;res@gsnLeftString = "m/s"
 resR@tiYAxisString = "Surface wind speed (m/s)"
 ;res@gsnCenterStringOrthogonalPosF = 0.05

 ; GSM GPV
 lfiles = (/gpvdir+"/np_fcst_phy_"+init+".nc", \ ; climatological SST
           gpvdir+"_est/np_fcst_phy_"+init+".nc", \ ; estimated SST (operational)
           gpvdir+"_mgdsst/np_fcst_phy_"+init+".nc"/) ; MGDSST
 rfiles = (/gpvdir+"/np_fcst_asia_"+init+".nc", \ ; climatological SST
           gpvdir+"_est/np_fcst_asia_"+init+".nc", \ ; estimated SST (operational)
           gpvdir+"_mgdsst/np_fcst_asia_"+init+".nc"/) ; MGDSST
 ltxts  = (/"flsh_avr_core_"+init+".txt", \
           "flsh_avr_core_"+init+"_est.txt", \
           "flsh_avr_core_"+init+"_mgdsst.txt"/)
 rtxts  = (/"twspd_avr_core_"+init+".txt", \
           "twspd_avr_core_"+init+"_est.txt", \
           "twspd_avr_core_"+init+"_mgdsst.txt"/)
 colors   = (/"red","blue","green"/)
 do j = 0, dimsizes(lfiles)-1
  fl = addfile(lfiles(j), "r")
  ;fr = addfile(rfiles(j), "r")
  tl = fl->time
  ;tr = fr->time
  lat = fl->lat
  r = new(dimsizes(lat), typeof(lat))
  deg2r(lat, r)
  print(sprintf("%.1f",r))
  r@units = "km"
  r = abs(r - rmax)
  imax = minind(r)
  print(imax)
  print(r(imax))
  ldata = new(ndate, float)
  ;rdata = new(ndate, float)
  ;data(0) = default_fillvalue("float")
  ;time(0) = 0.0
  ;do t = 0, ndate-2
  ;  var = f->FLLH(t,:imax,:)
  ;  data(t+1) = avg(var)
  ;  time(t+1) = (t+1)*6.*60.*60.
  ;end do
  do t = 0, ndate-1
    if (t.gt.0) then
      var = fl->FLSH(t-1,1:imax+1,:)
      ldata(t) = avg(var)
    ;  ind_uv = ind(tr.eq.tl(t-1))
    else
      ldata(t) = default_fillvalue("float")
    ;  ind_uv = 0
    end if
    ;ind_uv = t
    ;u = fr->U10m(ind_uv,1:imax+1,:)
    ;;v = fr->V10m(ind_uv,1:imax+1,:)
    ;;rdata(t) = avg(sqrt(u^2+v^2))
    ;rdata(t) = avg(sqrt(u^2))
  end do
  ldata!0 = "time"
  ldata&time = time
  printMinMax(ldata, 0)
  asciiwrite(ltxts(j), ldata)
  ;rdata!0 = "time"
  ;rdata&time = time
  ;printMinMax(rdata, 0)
  ;asciiwrite(rtxts(j), rdata)
  resR@xyLineColor = colors(j)
  resL@xyLineColor = colors(j)
  ; res@tmLabelAutoStride = True
  ;if (j.eq.0) then
  ;  plot = gsn_csm_xy2(wks, time, ldata, rdata, resL, resR)
  ;else
  ;  plot2 = gsn_csm_xy2(wks, time, ldata, rdata, resL, resR)
  ;  overlay(plot, plot2)
  ;  delete(plot2)
  ;end if
  ;data = data / radar
  ;print(data)
  delete(ldata)
  ;delete(rdata)
 end do
 ; Legends
 gres = True
 ;gres@YPosPercent = 95. ;upper
 gres@YPosPercent = 15. ;lower
 gres@XPosPercent = 5. ;left
 lres = True
 lres@lgLineColors = colors
 lres@lgLineThicknesses = 2.5
 lres@LineLengthPercent = 9.
 tres = True
 tres@lgLabels = (/"clim", "est", "mgdsst"/)
 tres@lgLineLabelFontHeightF = 0.1
 ;plotl = simple_legend(wks, plot, gres, lres, tres)
 ;draw(plot)
 ;frame(wks)
end