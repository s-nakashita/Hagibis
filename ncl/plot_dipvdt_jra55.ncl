load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
  dirname = "/Volumes/dandelion/grib/dias/data/jra55/Hist/Daily/anl_isentrop125/201910"
  tunits = "hours since 2019-10-05 00:00:00"
  time = ispan(0, 210, 6)
  time@units = tunits
  ntime= dimsizes(time)
  level = 340 ;K
  navg = 4
  iplot = 1
  lfilter = False
  fil = ""
  ;; d(PV)/dt
  do it = iplot,iplot+navg-1 ;ntime-2
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
  tmp = a->pVOR_GDS0_THEL({level},::-1,:)
  tmp = tmp * 1.0e6
  tmp@units = "PVU"
  if (it.eq.iplot) then
  pv = tmp
  idate = pv@initial_time
  else
  pv = pv + tmp
  edate = tmp@initial_time
  end if
  delete(a)
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
  if (it.eq.iplot) then
  ;uavg = u 
  ;vavg = v
  gpvx = tmp
  gpvy = tmp
  gradsf(tmp,gpvx,gpvy)
  ;gradlist = grad_latlon_cfd(tmp, lat, lon, True, False)
  ;gpvy = gradlist[0]
  ;gpvx = gradlist[1]
  gpvx = u*gpvx
  gpvy = v*gpvy
  gpvx@long_name = "Longitudinal advection of potential vorticity"
  gpvy@long_name = "Latitudinal advection of potential vorticity"
  gpvx@units = "PVU/m"
  gpvy@units = "PVU/m"
  else
  ;uavg = uavg + u 
  ;vavg = vavg + v
  tmpgpvx = tmp
  tmpgpvy = tmp
  gradsf(tmp,tmpgpvx,tmpgpvy)
  ;gradlist = grad_latlon_cfd(tmp, lat, lon, True, False)
  ;tmpgpvy = gradlist[0]
  ;tmpgpvx = gradlist[1]
  gpvx = gpvx + u*tmpgpvx
  gpvy = gpvy + v*tmpgpvy
  end if
  ;;;before
  date = cd_string(time(it-1),"%Y%N%D%H")
  vname = "pvor"
  fname = dirname + "/anl_isentrop125_"+vname+"."+date
  a = addfile(fname,"r")
  tmp1   = a->pVOR_GDS0_THEL({level},::-1,:)
  tmp1 = tmp1 * 1.0e6
  tmp1@units = "PVU/s"
  delete(a)
  ;;;after
  date = cd_string(time(it+1),"%Y%N%D%H")
  vname = "pvor"
  fname = dirname + "/anl_isentrop125_"+vname+"."+date
  a = addfile(fname,"r")
  tmp2   = a->pVOR_GDS0_THEL({level},::-1,:)
  tmp2 = tmp2 * 1.0e6
  tmp2 = (tmp2 - tmp1) / (12.*3600.)
  if (it.eq.iplot) then
  dpv = tmp2
  else
  dpv = dpv + tmp2
  end if
  delete(a)
  end do
  ;; advection term
  ;uavg = uavg / navg
  ;vavg = vavg / navg
  gpvx = gpvx / navg
  gpvy = gpvy / navg
  pv = pv / navg
  dpv = dpv / navg
;*************************************************
; hoskins filtering
;*************************************************
  if (lfilter) then
    fil = "fil_"
    ndim = dimsizes(pv)
    nlat = ndim(0)
    nlon = ndim(1)

    N = nlat
    if (nlon%2 .eq. 0)then
      N = min((/ nlat, (nlon+2)/2 /))
    else
      N = min((/ nlat, (nlon+1)/2 /))
    end if

    T = 42
    ;NT = 63
    a = new((/ nlat, N /),float)
    b = new((/ nlat, N /),float)
    ;; pv
    tmp = pv
    shaec(tmp, a, b)
    tri_trunc(a, b, T)
    ;exp_tapersh(a, b, NT, 2)
    shsec(a, b, tmp)
    pv = (/tmp/)
    ;; dpv
    tmp = dpv
    shaec(tmp, a, b)
    tri_trunc(a, b, T)
    ;exp_tapersh(a, b, NT, 2)
    shsec(a, b, tmp)
    dpv = (/tmp/)
    ;;; u
    ;tmp = uavg
    ;shaec(tmp, a, b)
    ;tri_trunc(a, b, T)
    ;;exp_tapersh(a, b, NT, 2)
    ;shsec(a, b, tmp)
    ;uavg = (/tmp/)
    ;;; v
    ;tmp = vavg
    ;shaec(tmp, a, b)
    ;tri_trunc(a, b, T)
    ;;exp_tapersh(a, b, NT, 2)
    ;shsec(a, b, tmp)
    ;vavg = (/tmp/)
    ;; gpvx
    tmp = gpvx
    shaec(tmp, a, b)
    tri_trunc(a, b, T)
    ;exp_tapersh(a, b, NT, 2)
    shsec(a, b, tmp)
    gpvx = (/tmp/)
    ;; gpvy
    tmp = gpvy
    shaec(tmp, a, b)
    tri_trunc(a, b, T)
    ;exp_tapersh(a, b, NT, 2)
    shsec(a, b, tmp)
    gpvy = (/tmp/)
    delete(a)
    delete(b)
  end if
  ;gpvy = pv
  ;gpvx = pv
  ;gradsf(pv, gpvx, gpvy)
  ;;gradlist = grad_latlon_cfd(pv, pv&g0_lat_1, pv&g0_lon_2, True, False)
  ;;gpvy = gradlist[0]
  ;;gpvx = gradlist[1]
  ;gpvx = uavg*gpvx
  ;gpvy = vavg*gpvy
  ;gpvx@long_name = "Longitudinal advection of potential vorticity"
  ;gpvy@long_name = "Latitudinal advection of potential vorticity"
  ;gpvx@units = "PVU/s"
  ;gpvy@units = "PVU/s"
  adv = pv
  adv = gpvx + gpvy
  adv@long_name = "Horizontal advection of potential vorticity"
  adv@units = "PVU/s"
  printVarSummary(pv)
  printMinMax(pv, 0)
  printVarSummary(dpv)
  printMinMax(dpv, 0)
  printVarSummary(gpvx)
  printMinMax(gpvx, 0)
  printVarSummary(gpvy)
  printMinMax(gpvy, 0)
  printVarSummary(adv)
  printMinMax(adv,0)
  cmap = read_colormap_file("GMT_polar")
  figname = "dIPVdt_np_"+sprinti("%0.3d",level)+"K_"+fil
  if (navg.eq.1) then
  figname = figname+sprinti("%d",iplot)
  else
  figname = figname+"navg"+sprinti("%d",iplot)+"-"+sprinti("%d",iplot+navg-1)
  end if
  wks = gsn_open_wks("png",figname)
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
  res@gsnPolar = "NH"
  ;res@mpMinLonF            =  110.               ; select a subregion
  ;res@mpMaxLonF            =  160.
  res@mpMinLatF            =  40. 
  ;res@mpMaxLatF            =  60.

  res@lbLabelBarOn         = True
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -9.               ; set min contour level
  res@cnMaxLevelValF       = 9.                ; set max contour level
  res@cnLevelSpacingF      = 1.                ; set contour spacing
  
  res@gsnLeftString = "PV"
  ;plots(0) = gsn_csm_contour_map_ce(wks,pv,res)
  plots(0) = gsn_csm_contour_map_polar(wks,pv,res)

  ;res@lbLabelBarOn         = False
  ;res@cnLevelSelectionMode = "ManualLevels"
  ;res@cnMinLevelValF       = -25.0               ; set min contour level
  ;res@cnMaxLevelValF       = 25.0                ; set max contour level
  ;res@cnLevelSpacingF      = 5.0                 ; set contour spacing
  ;res@gsnRightString = "m/s"
  ;res@gsnLeftString = "u"
  ;;plots(1) = gsn_csm_contour_map_ce(wks,gpvx,res)
  ;plots(1) = gsn_csm_contour_map_polar(wks,uavg,res)
  ;res@gsnLeftString = "v"
  ;;plots(2) = gsn_csm_contour_map_ce(wks,gpvy,res)
  ;plots(2) = gsn_csm_contour_map_polar(wks,vavg,res)
  ;delete(res@gsnRightString)
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -20.0               ; set min contour level
  res@cnMaxLevelValF       = 20.0                ; set max contour level
  res@cnLevelSpacingF      = 2.0                 ; set contour spacing
  res@gsnRightString = "10^-5 PVU/s"
  res@gsnLeftString = "u*dPdx"
  gpvx = gpvx * 1.0e5
  ;plots(1) = gsn_csm_contour_map_ce(wks,gpvx,res)
  plots(1) = gsn_csm_contour_map_polar(wks,gpvx,res)
  res@gsnLeftString = "v*dPdy"
  gpvy = gpvy * 1.0e5
  ;plots(2) = gsn_csm_contour_map_ce(wks,gpvy,res)
  plots(2) = gsn_csm_contour_map_polar(wks,gpvy,res)
  delete(res@gsnRightString)
  res@gsnLeftString = "Local time derivative"
  dpv = dpv * 1.0e5
  ;plots(3) = gsn_csm_contour_map_ce(wks,dpv,res)
  plots(3) = gsn_csm_contour_map_polar(wks,dpv,res)
  res@gsnLeftString = "Horizontal advection"
  adv = adv * 1.0e5
  ;plots(4) = gsn_csm_contour_map_ce(wks,adv,res)
  plots(4) = gsn_csm_contour_map_polar(wks,adv,res)
  res@gsnLeftString = "Sum"
  dpv = dpv + adv
  ;plots(5) = gsn_csm_contour_map_ce(wks,dpv,res)
  plots(5) = gsn_csm_contour_map_polar(wks,dpv,res)

  resP	                   = True
  title = sprinti("%0.3d",level)+"K "+idate
  if (navg .gt. 1) then
    title = title+"-"+edate
  end if
  resP@gsnPanelMainString	 = title
  ;resP@gsnPanelLabelBar    = True
  gsn_panel(wks,plots,(/2,3/),resP)
  delete(plots)
end