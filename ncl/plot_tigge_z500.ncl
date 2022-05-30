;lon0 = 120
;lon1 = 150
;lat0 =  0
;lat1 =  80
;fh = 0
;yyyymmddhh = 2019101200
;dprcp = 50
;dev = "png"
;orig = ecmf
begin
  center = "jma"
  yyyy = str_get_cols(sprinti("%0.10d", yyyymmddhh), 0, 3)
  mm = str_get_cols(sprinti("%0.10d", yyyymmddhh), 4, 5)
  dd = str_get_cols(sprinti("%0.10d", yyyymmddhh), 6, 7)
  hh = str_get_cols(sprinti("%0.10d", yyyymmddhh), 8, 9)
  fm = addfile("/Volumes/dandelion/netcdf/tigge/"+yyyy+"/"+center+"/glb_"+yyyymmddhh+"_mean.nc", "r")
  fs = addfile("/Volumes/dandelion/netcdf/tigge/"+yyyy+"/"+center+"/glb_"+yyyymmddhh+"_spr.nc", "r")
;  date_slp = fs->Date
;  date_ap = fp->Date
;  ind_slp = ind(date_slp.eq.2019-10-12)
;  ind_ap = ind(date_ap.eq.2019-10-12)
;  f = addfile("../hourly/"+yyyy+"/"+mm+"/"+dd+"/jrap"+yyyymmddhh+".nc", "r")

  wks = gsn_open_wks(dev, "tigge_"+center+"_pl_"+yyyymmddhh+"_h"+fh)
;  gsn_define_colormap(wks, "precip2_17lev")
  lon0 = -30
  lon1 = 210
  lat0 =  -40
  lat1 =  80

  res1 = True
  res1@gsnDraw = False
  res1@gsnFrame = False
  res1@gsnAddCyclic = True
  res1@mpDataBaseVersion = "MediumRes"
;  res1@mpLimitMode = "LatLon"
  res1@mpGeophysicalLineThicknessF = 2
  res1@mpMinLonF = lon0
  res1@mpMaxLonF = lon1
  res1@mpCenterLonF = 0.5 * (lon0 + lon1)
  res1@mpMinLatF = lat0
  res1@mpMaxLatF = lat1
  res1@cnSmoothingOn = True
  res1@cnFillMode = "RasterFill"
;  fillcol = ispan(2, 19, 1)
;  fillcol(0) = -1
;  res1@cnFillColors = fillcol
  res1@cnRasterSmoothingOn = True
  res1@cnLevelSelectionMode = "ExplicitLevels"
  res1@cnLevels = ispan(1,10,1)
;  res1@cnMinlevelValF = 0.
;  res1@cnMaxlevelValF = 50.
;  res1@cnLevelSpacingF = 10.
  res1@cnFillOn = True
  res1@cnFillPalette = "WhiteBlueGreenYellowRed"
  res1@cnLinesOn = False
  res1@cnLineLabelsOn = False
;  res@mpDataBaseVersion = "Ncarg4_1"
  res1@tiMainString = str_upper(center)+" z500+spread init "+yyyymmddhh ;201910"+sprinti("%0.2d",id+d)
  res1@gsnLeftString = " "
;  res1@gsnLeftString = "sea level pressure hPa"
;  res1@gsnLeftString = "geo-potential height gpm"
  res1@gsnRightString = "lead time (hour)="+fh
  res1@lbOrientation = "vertical"
  res1@vpXF = 0.15
  res1@vpWidthF = 0.65


  res2 = True
  res2@gsnDraw = False
  res2@gsnFrame = False
  res2@gsnAddCyclic = True
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
;  res2@cnLevels = (/948, 952, 956, 960, 964, 968, 972, 976, 980, 984, 988, 992, 996, 1000, 1004, 1008, 1012, 1016, 1020/)
;   res2@cnLevels = (/5200, 5250, 5300, 5350, 5400, 5450, 5500, 5550, 5600, 5650, 5700, 5750, 5800, 5850, 5900, 5950, 6000/)
  res2@cnLevels = ispan(4800,5840,80)
;  else
;    res1@cnLevels = dprcp * ispan(1, 17, 1)
;  end if
  res2@cnLineLabelsOn = False
  
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
  res2@gsnLeftString = ""
  res2@gsnRightString = ""
;  res@lbOrientation = "vertical"
;  res@vpXF = 0.15
;  res@vpWidthF = 0.65


;  n = 10
;  lon_in = fm->longitude
;  lon_out = lon_in(n/2::n)
;  lat_in = fm->latitude
;  lat_out = lat_in(n/2::n)
  it = fh / 6
  print(it)
  z500 = fm->gh(it,{500},:,:)
  spr = fs->gh(it,{500},:,:)
;  printVarSummary(spr)
;  deg2rad = acos(-1) / 180.
;  slp = area_hi2lores_Wrap(lon_in, lat_in, slp_in, False, cos(lat_in * deg2rad), lon_out, lat_out, False) 
;  spr = area_hi2lores_Wrap(lon_in, lat_in, spr_in, False, cos(lat_in * deg2rad), lon_out, lat_out, False) 
;  slp = slp*0.01
;  spr = spr*0.001
  plot = gsn_csm_contour_map(wks,spr,res1)
  plot2 = gsn_csm_contour(wks,z500,res2)
  overlay(plot, plot2)
  draw(plot)
  frame(wks)
;  plot = gsn_csm_contour_map_overlay(wks, spr, z500, res1, res2)	
;  plot = gsn_csm_contour_map(wks, slp, res1)

end
