;lon0 = 120
;lon1 = 150
;lat0 =  0
;lat1 =  80
;yyyymmdd = 20191011
;d = 0
;yyyymmddhh = 2019101200
;dprcp = 50
;dev = "png"
;orig = ecmf
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
  yyyymmddhh = 2019100912
  yyyy = str_get_cols(sprinti("%0.10d", yyyymmddhh), 0, 3)
  mmddhh = str_get_cols(sprinti("%0.10d", yyyymmddhh), 4, 9)
  mm = str_get_cols(sprinti("%0.10d", yyyymmddhh), 4, 5)
  dd = str_get_cols(sprinti("%0.10d", yyyymmddhh), 6, 7)
  hh = str_get_cols(sprinti("%0.10d", yyyymmddhh), 8, 9)
  id = tointeger(dd)
  fm = addfile("/Users/nakashita/netcdf/tigge/"+yyyy+"/jma/"+mmddhh+"_mean.nc", "r")
  fs = addfile("/Users/nakashita/netcdf/tigge/"+yyyy+"/jma/"+mmddhh+"_spr.nc", "r")
  time  = fm->time
  date  = cd_string(time, "%Y%N%D%H")
;  date_slp = fs->Date
;  date_ap = fp->Date
;  ind_slp = ind(date_slp.eq.2019-10-12)
;  ind_ap = ind(date_ap.eq.2019-10-12)
;  f = addfile("../hourly/"+yyyy+"/"+mm+"/"+dd+"/jrap"+yyyymmddhh+".nc", "r")

  lon0 = 120.
  lon1 = 160.
  lat0 =  10.
  lat1 =  60.

  wks = gsn_open_wks("png", "tigge_jma_slp"+yyyymmddhh)
  plot = new(8,graphic)
;  gsn_define_colormap(wks, "precip2_17lev")
  vdate = yyyymmddhh
  do d=0,7
    print(vdate)
    index   = ind(date.eq.vdate)
    print(index)
    res1 = True
    res1@gsnDraw	      = False
    res1@gsnFrame	      = False
    res=False

    res1@gsnAddCyclic = False
    res1@mpDataBaseVersion = "MediumRes"
    res1@mpLimitMode = "LatLon"
    res1@mpGeophysicalLineThicknessF = 2
    res1@mpMinLonF = lon0
    res1@mpMaxLonF = lon1
    res1@mpCenterLonF = 0.5 * (lon0 + lon1)
    res1@mpMinLatF = lat0
    res1@mpMaxLatF = lat1
    res1@cnSmoothingOn = True
;  res1@cnFillMode = "RasterFill"
;  fillcol = ispan(2, 19, 1)
;  fillcol(0) = -1
;  res1@cnFillColors = fillcol
    res1@cnRasterSmoothingOn = True
    res1@cnLevelSelectionMode = "ExplicitLevels"
    res1@cnLevels = (/2., 4., 6., 8., 10., 12., 14., 16., 18., 20./)
;  res1@cnMinlevelValF = 0.
;  res1@cnMaxlevelValF = 50.
;  res1@cnLevelSpacingF = 10.
    res1@cnFillOn = True
    res1@cnFillPalette = "WhiteBlueGreenYellowRed"
    res1@cnLinesOn = False
    res1@cnLineLabelsOn = False
;  res@mpDataBaseVersion = "Ncarg4_1"
    res1@tiMainString = "tigge_jma slp+spread "+sprinti("%0.10d",vdate)
    res1@gsnLeftString = "hPa"
;  res1@gsnLeftString = "sea level pressure hPa"
;  res1@gsnLeftString = "geo-potential height gpm"
    res1@gsnRightString = "initial date="+yyyymmddhh
    res1@lbOrientation = "vertical"
    res1@vpXF = 0.15
    res1@vpWidthF = 0.65


    res2 = True
    res2@gsnDraw	      = False
    res2@gsnFrame	      = False
;  res=False

    res2@gsnAddCyclic = False
;  res2@mpMinLonF = lon0
;  res2@mpMaxLonF = lon1
;  res2@mpCenterLonF = 0.5 * (lon0 + lon1)
;  res2@mpMinLatF = lat0
;  res2@mpMaxLatF = lat1
    res2@cnSmoothingOn = True
;  res2@cnFillOn = False
;  res@cnFillMode = "RasterFill"
;  fillcol = ispan(2, 19, 1)
;  fillcol(0) = -1
;  res@cnFillColors = fillcol
;  res@cnRasterSmoothingOn = True
;  res2@cnLinesOn = True
    res2@cnLevelSelectionMode = "ExplicitLevels"
;  res2@cnMinLevelValF = 960.
;  res2@cnMaxLevelValF = 1020.
;  res2@cnLevelSpacingF = 4.  
;  if dprcp .eq. 0 then
    res2@cnLevels = (/948, 952, 956, 960, 964, 968, 972, 976, 980, 984, 988, 992, 996, 1000, 1004, 1008, 1012, 1016, 1020/)
;    res2@cnLevels = (/5200, 5250, 5300, 5350, 5400, 5450, 5500, 5550, 5600, 5650, 5700, 5750, 5800, 5850, 5900, 5950, 6000/)
;  else
;    res1@cnLevels = dprcp * ispan(1, 17, 1)
;  end if
  
    res2@cnLineColor = "blue"
;  if dprcp .eq. 0 then
;    res@cnLevels = (/10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300/)
;  else
;    res@cnLevels = dprcp * ispan(1, 17, 1)
;  end if
;  res@mpDataBaseVersion = "Ncarg4_1"
;  res2@mpDataBaseVersion = "RANGS_GSHHS"
;  res2@mpGeophysicalLineThicknessF = 2
;  res@tiMainString = "GSM around Japan"
;  res@gsnLeftString = "daily prcp mm"
;  res2@gsnLeftString = "sea level pressure Pa"
;  res@gsnRightString = yyyymmdd
;  res@lbOrientation = "vertical"
;  res@vpXF = 0.15
;  res@vpWidthF = 0.65


;  n = 10
;  lon_in = fm->longitude
;  lon_out = lon_in(n/2::n)
;  lat_in = fm->latitude
;  lat_out = lat_in(n/2::n)
    slp = fm->PRES_meansealevel(index,:,:)
    spr = fs->PRES_meansealevel(index,:,:)
;  printVarSummary(spr)
;  deg2rad = acos(-1) / 180.
;  slp = area_hi2lores_Wrap(lon_in, lat_in, slp_in, False, cos(lat_in * deg2rad), lon_out, lat_out, False) 
;  spr = area_hi2lores_Wrap(lon_in, lat_in, spr_in, False, cos(lat_in * deg2rad), lon_out, lat_out, False) 
    slp = slp*0.01
    spr = spr*0.01
    plot(d) = gsn_csm_contour_map_overlay(wks, spr, slp, res1, res2)	
;  plot = gsn_csm_contour_map(wks, slp, res1)
    year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
    month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
    day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
    hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
    sec = cd_inv_calendar(year,month,day,hour,0,0, time@units, 0)
    newsec = sec
    newsec = sec + 12*60*60
    vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
  end do
  resP	                   = True
  resP@gsnPanelMainString	   = "SLP+spread JMA ensemble forecast from "+sprinti("%0.10d",yyyymmddhh)
  gsn_panel(wks,plot,(/2,4/),resP)
end
