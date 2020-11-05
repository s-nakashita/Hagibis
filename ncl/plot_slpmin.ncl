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
 
  lon0 = 139.0
  lat0 = 22.5
  lonc = 139.05
  latc = 22.6
  
  wks = gsn_open_wks("png", "tigge_jma_slp"+yyyymmddhh)
;  plot = new(8,graphic)
;  gsn_define_colormap(wks, "precip2_17lev")
  vdate = 2019101000
;  do d=2, 2
  print(vdate)
  index   = ind(date.eq.vdate)
  print(index)
  res1 = True
  res1@gsnDraw	      = False
  res1@gsnFrame	      = False
  res1@gsnAddCyclic = False
  res1@mpDataBaseVersion = "MediumRes"
  res1@mpLimitMode = "LatLon"
  res1@mpGeophysicalLineThicknessF = 2
  res1@mpMinLonF = lon0 - 1.
  res1@mpMaxLonF = lon0 + 1.
  res1@mpCenterLonF = lon0
  res1@mpMinLatF = lat0 - 1.
  res1@mpMaxLatF = lat0 + 1.
  res1@cnSmoothingOn = True
;  res@mpDataBaseVersion = "Ncarg4_1"
  res1@tiMainString = "tigge_jma slp-min "+sprinti("%0.10d",vdate)
  res1@gsnLeftString = "hPa"
  res1@lbOrientation = "vertical"
  res1@vpXF = 0.15
  res1@vpWidthF = 0.65

  ;res1@cnSmoothingOn = True
  ;res1@cnRasterSmoothingOn = True
  res1@cnLevelSelectionMode = "ExplicitLevels"
  res1@cnLevels = (/960, 964, 968, 972, 976, 980, 984, 988, 992, 996, 1000, 1004/)
  
;    res1@cnLineColor = "blue"

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
  plot  = gsn_csm_contour_map(wks, slp, res1)	

  pmres = True
  pmres@gsMarkerColor = "blue"
  pmres@gsFillColor = "blue"
  pmres@gsMarkerSizeF = 0.03
  pmres@gsLineThicknessF = 5.
  pmres@gsMarkerIndex = 12
  str = unique_string("poly")
  x = lonc
  y = latc
  plot@$str$ = gsn_add_polymarker(wks, plot, x, y, pmres)

  pmres@gsMarkerColor = "red"
  pmres@gsFillColor = "red"
  pmres@gsMarkerSizeF = 0.02
  pmres@gsLineThicknessF = 3.
  pmres@gsMarkerIndex = 5
  do i = -1, 1
    do j = -1, 1
      str = unique_string("poly")
      x = lon0 + i*0.5
      y = lat0 + j*0.5
      plot@$str$ = gsn_add_polymarker(wks, plot, x, y, pmres)
    end do
  end do
  draw(plot)
  frame(wks)
end
