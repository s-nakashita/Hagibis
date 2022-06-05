load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
 rdrdir = "/Volumes/dandelion/netcdf/JMAGPV/hourly/2019/10"
 gpvdir = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"
 init = "2019100912"

 region = "region1"
 wks = gsn_open_wks("png", "prcp_"+region+"_"+init)
 ;; hakone
 ;lonh = 139
 ;lath = 35
 ;lonmin = 138.875
 ;lonmax = 139.125
 ;latmin = 34.875
 ;latmax = 35.125
 if (region .eq. "region1") then
 ;;region1 (Tokai)
 lonmin = 136.
 lonmax = 138.
 latmin = 34.
 latmax = 36.
 ;;region2
 ;lonmin = 137.5
 ;lonmax = 141.5
 ;latmin = 34.5
 ;latmax = 38.5
 else if (region .eq. "KN") then
 ;;Kanto
 lonmin = 138.
 lonmax = 141.
 latmin = 34.5
 latmax = 36.5
 else if (region .eq. "ST") then
 ;; Southern Tohoku
 lonmin = 139.
 lonmax = 141.5
 latmin = 36.5
 latmax = 38.5
 end if
 end if
 end if

 res = True
 res@gsnDraw = True
 res@gsnFrame = False

 ;accumulated value plotting resources
 xyres = res
 xyres@trYMinF = 0.0
 xyres@trYMaxF = 300.0
 xyres@tmXBOn  = False
 xyres@tmXTOn  = False
 xyres@tmYLOn  = False
 xyres@tmYROn  = True
 xyres@tmYRLabelsOn  = True

 xyres@tiYAxisString = "Accumulated precipitation (mm)"
 xyres@tiYAxisSide = "Right"
 xyres@tiYAxisAngleF = 270

 ;intensity plotting resources
 bres = res 
 bres@gsnXYBarChart = True
 bres@tmYROn = False
 bres@trYMinF = 0.0
 bres@trYMaxF = 30.0
 bres@tiYAxisString = "Hourly precipitation (mm)"
 bres@tiYAxisSide   = "Left"
 bres@tiYAxisAngleF = 90
 ;res@gsnCenterString = "3h precipitation at "+sprinti("%d", lonh)+"E, "+sprinti("%d", lath)+"N"
 bres@gsnCenterString = sprintf("%.1f", lonmin)+"E-"+sprintf("%.1f", lonmax)+\
 "E, "+sprintf("%.1f", latmin)+"N-"+sprintf("%.1f", latmax)+"N"
 bres@gsnCenterStringOrthogonalPosF = 0.05

 ndate = 47
 ;ndate = 15
 data = new((/ndate,4/), float) ; precipitation intensity (mm/hr)
 adata = new((/ndate,4/), float) ; accumulated precipitation (mm)
 time = new(ndate, float)
 time@units = "seconds since 2019-10-11 00:00:00"
 tmp = 0.
 delta = 3600.;*3
 data(:,:) = default_fillvalue("float")
 adata(:,:) = default_fillvalue("float")
 ; Radar
 do t = 0, ndate-1
   ;tmp = 0.
   ;do h = 1, 3
   h = 1
     if (t+h .lt. 24) then
        hh = sprinti("%0.2d", t+h)
        fname = rdrdir + "/11/jrap20191011"+hh+".nc"
     else
        hh = sprinti("%0.2d", t+h-24)
        fname = rdrdir + "/12/jrap20191012"+hh+".nc"
     end if
     f = addfile(fname, "r")
     prcp_in = f->var0_1_201_surface(0,{latmin:latmax},{lonmin:lonmax})
     ;printMinMax(prcp_in, 0)
     tmp = tmp + avg(prcp_in)
   ;end do
   if (t.eq.0) then
     data(t,0) = tmp ;/ 3.0
   else
     data(t,0) = (tmp - adata(t-1,0));/3.0
   end if
   adata(t,0) = tmp
   time(t) = (t+1)*delta
 end do
 data!0 = "time"
 data&time = time
 printMinMax(data, 0)
 adata!0 = "time"
 adata&time = time
 printMinMax(adata, 0)
 
 bres@tmXBMode = "Explicit"
 bres@tmXBValues = time(5::6)
 bres@tmXBLabels = cd_string(time(5::6), "%D-%HZ")
 bres@tmXBMinorOn = True
 bres@tmXBMinorValues = time
 bres@trXMinF = time(0) - 0.5*delta
 bres@trXMaxF = time(ndate-1) + 0.5*delta
 xyres@trXMinF = time(0) - 0.5*delta
 xyres@trXMaxF = time(ndate-1) + 0.5*delta
; res@tmLabelAutoStride = True

 ; GSM GPV
 files = (/gpvdir+"/fcst_surf_asia_"+init+".nc", \ ; climatological SST
           gpvdir+"_est/fcst_surf_asia_"+init+".nc", \ ; estimated SST (operational)
           gpvdir+"_mgdsst/fcst_surf_asia_"+init+".nc"/) ; MGDSST
 colors   = (/"red","blue","green"/)
 dum = new(3, graphic)
 do j = 0, dimsizes(files)-1
  fname = files(j)
  f = addfile(fname, "r")
  ftime = f->time
  nftime = dimsizes(ftime)
  tunits = ftime@units
  tfirst = cd_inv_calendar( 2019, 10, 11, 1, 0, 0, tunits, 0)
  ind_first = ind(f->time.eq.tfirst)
  do t = 0, ndate-1
    if (ind_first+t .gt. nftime-1) then
      break
    end if
    if (t .lt. 24) then
      t0h = cd_inv_calendar( 2019, 10, 11, t, 0, 0, tunits, 0)
      if (t .eq. 23) then
        t3h = cd_inv_calendar( 2019, 10, 12, 0, 0, 0, tunits, 0)
      else
        t3h = cd_inv_calendar( 2019, 10, 11, t+1, 0, 0, tunits, 0)
      end if
    else
      t0h = cd_inv_calendar( 2019, 10, 12, t-24, 0, 0, tunits, 0)
      if (t .eq. 48) then
        t3h = cd_inv_calendar( 2019, 10, 13, 0, 0, 0, tunits, 0)
      else
        t3h = cd_inv_calendar( 2019, 10, 12, t-23, 0, 0, tunits, 0)
      end if
    end if
    ;prcp1 = f->RAIN({t0h},{lath},{lonh})
    ;prcp2 = f->RAIN({t3h},{lath},{lonh})
    ;data(t) = (/prcp2 - prcp1/)
    prcp1 = f->RAIN({t0h},{latmin:latmax},{lonmin:lonmax})
    prcp2 = f->RAIN({t3h},{latmin:latmax},{lonmin:lonmax})
    data(t,j+1) = avg(prcp2-prcp1);/3.0
    adata(t,j+1) = avg(prcp2)
  end do
 printMinMax(data, 0)
 printMinMax(adata, 0)
 end do
 ;;Draw Bar plot
 ;;Radar
 bres@gsnXYBarChartBarWidth = 0.2*delta
 bres@gsnXYBarChartColors = "gray"
 base = time
 base = time - conform(time, bres@gsnXYBarChartBarWidth*1.5, -1)
 print(base)
 bar_plot = gsn_csm_xy(wks, base,  data(:,0), bres)
 base = base + conform(base, bres@gsnXYBarChartBarWidth, -1)
 ;;GSM GPV
 ;bres@gsnXYBarChartBarWidth = 0.1*delta
 bres@gsnXYBarChartFillOpacityF = 0.5
 do j = 1, dimsizes(files)
  bres@gsnXYBarChartColors = colors(j-1)
  bar_plot = gsn_csm_xy(wks, base,  data(:,j), bres)
  base = base + conform(base, bres@gsnXYBarChartBarWidth, -1)
 end do

 ;;Draw XY plot
 ;;Radar
 xyres@xyLineThicknessF = 6.0
 xyres@xyLineColor = "black"
 getvalues bar_plot
   "vpXF"   : vpx
   "vpYF"   : vpy 
   "vpWidthF" : vpw 
   "vpHeightF" : vph 
 end getvalues
 xyres@vpXF = vpx
 xyres@vpYF = vpy
 xyres@vpWidthF = vpw
 xyres@vpHeightF = vph
 xy_plot  = gsn_csm_xy(wks, time, adata(:,0), xyres)
 ;;GSM GPV
 xyres@xyLineThicknessF = 4.0
 do j = 1, dimsizes(files)
  xyres@xyLineColor = colors(j-1)
  xy_plot  = gsn_csm_xy(wks, time, adata(:,j), xyres)
  ratio = data(:,j) 
  ratio = ratio / data(:,0)
  print(ratio)
 end do 
 ;; save data
 header = "    time (intensity) radar  clim   est   mgd (accumulate) radar  clim   est   mgd"
 wfile = "prcp_"+region+"_"+init+".txt"
 system("rm -f "+wfile)
 wdata = new((/ndate,9/), float)
 wdata(:,0) = time
 wdata(:,1:4) = data
 wdata(:,5:) = adata
 option = True
 option@row = False
 option@tspace = 0
 option@fout = wfile
 option@title = header
 write_matrix(wdata, "f8.0, 8e16.7", option)
 delete(data)
 delete(adata)
;  delete(bar_plot2)
;  delete(xy_plot2)
 ; Legends
 gres = True
 gres@YPosPercent = 95.
 gres@XPosPercent = 5.
 lres = True
 lres@lgLineColors = colors
 lres@lgLineThicknesses = 2.5
 lres@LineLengthPercent = 9.
 tres = True
 tres@lgLabels = (/"clim", "est", "mgdsst"/)
 tres@lgLineLabelFontHeightF = 0.1
 plotl = simple_legend(wks, bar_plot, gres, lres, tres)

 anno_id = gsn_add_annotation(bar_plot, xy_plot, False)
 pres = True
 pres@gsnMaximize = True
 maximize_output(wks,pres)
 ;draw(plot)
 ;frame(wks)
end