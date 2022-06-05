load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
  dirname = "/Volumes/dandelion/grib/dias/data/jra55/Hist/Daily/anl_isentrop125/201910"
  tunits = "hours since 2019-10-05 00:00:00"
  time = ispan(0, 210, 6)
  time@units = tunits
  ntime= dimsizes(time)
  level = 320 ;K
  ;; d(PV)/dt
  do it = 1, ntime-2
  date = cd_string(time(it),"%Y%N%D%H")
  print(date)
  vname = "pvor"
  fname = dirname + "/anl_isentrop125_"+vname+"."+date
  a = addfile(fname,"r")
  ;dNames = getfiledimnames(a)
  ;print(dNames)
  ;vNames = getfilevarnames(a)
  ;print(vNames)
  lat   = a->g0_lat_1(::-1)
  lon   = a->g0_lon_2
  pv    = a->pVOR_GDS0_THEL({level},::-1,:) ;S=>N
  pv = pv * 1.0e6
  pv@units = "PVU"
  printVarSummary(pv)
  printMinMax(pv, 0)
  delete(a)
  ;;;before
  date = cd_string(time(it-1),"%Y%N%D%H")
  fname = dirname + "/anl_isentrop125_"+vname+"."+date
  a = addfile(fname,"r")
  dpv   = a->pVOR_GDS0_THEL({level},::-1,:)
  dpv = dpv * 1.0e6
  dpv@units = "PVU/s"
  delete(a)
  ;;;after
  date = cd_string(time(it+1),"%Y%N%D%H")
  fname = dirname + "/anl_isentrop125_"+vname+"."+date
  a = addfile(fname,"r")
  tmp   = a->pVOR_GDS0_THEL({level},::-1,:)
  tmp = tmp * 1.0e6
  dpv = (tmp - dpv) / (12.*3600.)
  printVarSummary(dpv)
  printMinMax(dpv, 0)
  delete(a)
  ;; advection term
  
  date = cd_string(time(it),"%Y%N%D%H")
  vname = "ugrd"
  fname = dirname + "/anl_isentrop125_"+vname+"."+date
  a = addfile(fname,"r")
  u = a->UGRD_GDS0_THEL({level},::-1,:) ;S=>N
  delete(a)
  vname = "vgrd"
  fname = dirname + "/anl_isentrop125_"+vname+"."+date
  a = addfile(fname,"r")
  v = a->VGRD_GDS0_THEL({level},::-1,:) ;S=>N
  delete(a)
  gpvx = pv
  gpvy = pv
  gradlist = grad_latlon_cfd(pv, lat, lon, True, False)
  gpvy = gradlist[0]
  gpvx = gradlist[1]
  gpvx = u*gpvx
  gpvy = v*gpvy
  gpvx@long_name = "Longitudinal advection of potential vorticity"
  gpvy@long_name = "Latitudinal advection of potential vorticity"
  gpvx@units = "PVU/m"
  gpvy@units = "PVU/m"
  printVarSummary(gpvx)
  printMinMax(gpvx, 0)
  printVarSummary(gpvy)
  printMinMax(gpvy, 0)
  adv = pv
  adv = gpvx + gpvy
  adv@long_name = "Horizontal advection of potential vorticity"
  adv@units = "PVU/s"
  printVarSummary(adv)
  printMinMax(adv,0)
  cmap = read_colormap_file("GMT_polar")
  iplot = it
  wks = gsn_open_wks("png","dIPVdt_"+sprinti("%0.3d",level)+"K_"+sprinti("%d",iplot))
  plots = new(6,graphic)
;---Create plot
  res                  = True               ; plot mods desired
  res@gsnDraw	      = False
  res@gsnFrame	      = False
  res@gsnAddCyclic	      = True
  res@cnFillOn             = True               ; turn on color for contours
  res@cnLinesOn            = True              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off contour line labels
  res@cnFillPalette        = cmap

  res@mpFillOn             = False 
  res@mpMinLonF            =  110.               ; select a subregion
  res@mpMaxLonF            =  160.
  res@mpMinLatF            =  10. 
  res@mpMaxLatF            =  60.

  res@lbLabelBarOn         = True
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -5.               ; set min contour level
  res@cnMaxLevelValF       = 5.                ; set max contour level
  res@cnLevelSpacingF      = .5               ; set contour spacing
  
  res@gsnLeftString = "PV"
  plots(0) = gsn_csm_contour_map_ce(wks,pv,res)

  ;res@lbLabelBarOn         = False
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -5.0e-4               ; set min contour level
  res@cnMaxLevelValF       = 5.0e-4                ; set max contour level
  res@cnLevelSpacingF      = 5.0e-5               ; set contour spacing
  res@gsnRightString = "PVU/s"
  res@gsnLeftString = "u*dPVdx"
  plots(1) = gsn_csm_contour_map_ce(wks,gpvx,res)
  res@gsnLeftString = "v*dPVdy"
  plots(2) = gsn_csm_contour_map_ce(wks,gpvy,res)
  delete(res@gsnRightString)
  res@gsnLeftString = "Local time derivative"
  plots(3) = gsn_csm_contour_map_ce(wks,dpv,res)
  res@gsnLeftString = "Horizontal advection"
  plots(4) = gsn_csm_contour_map_ce(wks,adv,res)
  res@gsnLeftString = "Sum"
  dpv = dpv + adv
  plots(5) = gsn_csm_contour_map_ce(wks,dpv,res)

  resP	                   = True
  resP@gsnPanelMainString	 = sprinti("%0.3d",level)+"K "+pv@initial_time
  ;resP@gsnPanelLabelBar    = True
  gsn_panel(wks,plots,(/2,3/),resP)
  delete(plots)
  end do
end