;----------------------------------------------------------------------
; stream_11.ncl 
;
; Concepts illustrated:
;   - Drawing streamlines colored by another field over a map
;----------------------------------------------------------------------
; Many thanks to Alan Brammer of the University at Albany who provided
; the initial version of this example which was modified to color
; the streamlines based on the field "p1". Alan also contributed the
; gsn_csm_streamline_scalar_xxxx functions which we added to V6.3.0.
;----------------------------------------------------------------------
; This example will only work with NCL V6.3.0 or later.
;----------------------------------------------------------------------
load "../rotate/librotate.ncl"
;
; rotation TC center to NP (irot=0) or NP to TC center (irot=1)
;
function rotate(u,v,lonc,latc,irot)
local lonin, latin, lonout, latout, \
      rlonc, rlatc, rlonin, rlatin, rlonout, rlatout, \
      xd, yd, zd, xdi, ydi, zdi, xdr, ydr, zdr, \
      ur, vr
begin
  ; set dimensions
  ndim = dimsizes(u)
  nrank = dimsizes(ndim)
  if(nrank.eq.3) then
    nlev = ndim(0)
    nlat = ndim(1)
    nlon = ndim(2)
    idlat = 1
    idlon = 2
  elseif(nrank.eq.2) then
    nlat = ndim(0)
    nlon = ndim(1)
    idlat = 0
    idlon = 1
  else
    print("Invalid rank of u")
    exit
  end if
  latin = u&lat
  lonin = u&lon
  print("latin")
  printMinMax(latin,0)
  print("lonin")
  printMinMax(lonin,0)
  print("======")
  ; prepare output coordinates
  rlatout = new(nlon*nlat, float)
  rlonout = new(nlon*nlat, float)
  rlonc = lonc * deg2rad
  rlatc = latc * deg2rad
  rlonin = lonin * deg2rad
  rlatin = latin * deg2rad
  if(irot.eq.0) then
    rotate_lonlat(rlonc, rlatc, rlonin, rlatin, rlonout, rlatout)
  elseif(irot.eq.1) then
    inv_rotate_lonlat(rlonc, rlatc, rlonin, rlatin, rlonout, rlatout)
  end if
  latout = rlatout * rad2deg
  lonout = rlonout * rad2deg
  print("latout")
  print(dimsizes(latout))
  printMinMax(latout,0)
  print("lonout")
  print(dimsizes(lonout))
  printMinMax(lonout,0)
  print("======")
  ; (u,v) => (xd,yd,zd)
  xd = u 
  yd = u 
  zd = u
  rlonmd = conform(xd, rlonin, idlon)
  rlatmd = conform(xd, rlatin, idlat)
  uv2xyzd(u, v, rlonmd, rlatmd, xd, yd, zd)
  if(any(ismissing(xd)).or.\
     any(ismissing(yd)).or.\
     any(ismissing(zd))) then
    print("Fatal: Missing values exist in xyzd")
    return
  end if
  print("xd")
  print(dimsizes(xd))
  printMinMax(xd,0)
  print("yd")
  print(dimsizes(yd))
  printMinMax(yd,0)
  print("zd")
  print(dimsizes(zd))
  printMinMax(zd,0)
  print("======")
  xdi = linint2_points(lonin, latin, xd, True, lonout, latout, 0)
  ydi = linint2_points(lonin, latin, yd, True, lonout, latout, 0)
  zdi = linint2_points(lonin, latin, zd, True, lonout, latout, 0)
  if(any(ismissing(xdi)).or.\
     any(ismissing(ydi)).or.\
     any(ismissing(zdi))) then
    print("Warning: Missing values occur after interpolation")
    datamsg = linmsg(xdi,0)
    xdi = (/datamsg/)
    datamsg = linmsg(ydi,0)
    ydi = (/datamsg/)
    datamsg = linmsg(zdi,0)
    zdi = (/datamsg/)
  end if
  print("xdi")
  print(dimsizes(xdi))
  printMinMax(xdi,0)
  print("ydi")
  print(dimsizes(ydi))
  printMinMax(ydi,0)
  print("zdi")
  print(dimsizes(zdi))
  printMinMax(zdi,0)
  print("======")
  ; rotate coordinates
  xdr = xdi
  ydr = ydi 
  zdr = zdi 
  rloncmd = conform(xdr, rlonc, -1)
  rlatcmd = conform(xdr, rlatc, -1)
  if(irot.eq.0) then
    tc2np(rloncmd, rlatcmd, \
      xdi, ydi, zdi, xdr, ydr, zdr)
  elseif(irot.eq.1) then
    np2tc(rloncmd, rlatcmd, \
      xdi, ydi, zdi, xdr, ydr, zdr)
  end if
  if(any(ismissing(xdr)).or.\
     any(ismissing(ydr)).or.\
     any(ismissing(zdr))) then
    print("Fatal: Missing values occur after rotation")
    return
  end if
  ; (xd,yd,zd) => (u,v)
  ur1d = xdr
  vr1d = xdr
  if(nrank.eq.3) then
    rlonrmd = conform(xdr,\
      ndtooned(conform_dims((/nlat,nlon/), rlonin, 1)),1)
  else
    rlonrmd = ndtooned(conform_dims((/nlat,nlon/), rlonin, 1))
  end if
  xyzd2uv(xdr,ydr,zdr,rlonrmd,ur1d,vr1d)
  ur = u
  vr = v
  if(nrank.eq.3) then
    do k=1,nlev
      ur(k,:,:) = (/reshape(ur1d(k,:),(/nlat,nlon/))/)
      vr(k,:,:) = (/reshape(vr1d(k,:),(/nlat,nlon/))/)
    end do
  elseif(nrank.eq.2) then
    ur = (/reshape(ur1d,(/nlat,nlon/))/)
    vr = (/reshape(vr1d,(/nlat,nlon/))/)
  end if
  return([/ ur, vr /])
end
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
  wtype = "pdf"
  init = "2019100912"
  prtb = "pn"
  it = 4
  plev0=850
  plev1=500
  if(prtb.eq."cntl") then
  ;ncout = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/fcst_p_lwr_"+init+".nc","r")
  ncout = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/fcst_p_"+init+".nc","r")
  track = readAsciiTable("../pytrack/track"+init+"_gsm_tl479_est.txt",7,"float",0)
  title = "CNTL"
  else
  ncout = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+"+prtb+"/fcst_p_"+init+".nc","r")
  track = readAsciiTable("../pytrack/track"+init+"_gsm_tl479_est+"+prtb+".txt",7,"float",0)
  if(prtb.eq."p") then
    title = "RIDGE+"
  else
    title = "RIDGE-"
  end if
  end if
  time = ncout->time
  ;pgrid = ncout->Z(it,{plev},::-1,:)
  if (plev0.ne.plev1) then
  clev=sprinti("%d",plev0)+"-"+sprinti("%d",plev1)
  u3d = ncout->U(it,{plev0:plev1},::-1,:)
  v3d = ncout->V(it,{plev0:plev1},::-1,:)
  plevs = u3d&level
  dp = dpres_plevel_Wrap(plevs, plev0, plev1, 0)
  wgt = dp / (plev0 - plev1)
  u3d = u3d * conform(u3d,wgt,0)
  v3d = v3d * conform(v3d,wgt,0)
  ugrid = u3d(0,:,:)
  ugrid = dim_sum_n(u3d, 0)
  vgrid = v3d(0,:,:)
  vgrid = dim_sum_n(v3d, 0)
  else
  plev=plev0
  clev=sprinti("%d",plev)
  ugrid = ncout->U(it,{plev},::-1,:)
  vgrid = ncout->V(it,{plev},::-1,:)
  end if
  ;rotate
  idy = toint(cd_string(time(it),"%D"))
  ihr = toint(cd_string(time(it),"%H"))
  npt = dimsizes(track(:,0))
  itc = 0
  dd = toint(track(:,2))
  hh = toint(track(:,3))
  do i = 0, npt-1
    if(dd(i).eq.idy.and.hh(i).eq.ihr) then
      break
    end if
  end do
  lonc = track(i,4)
  latc = track(i,5)
  print("center=("+lonc+","+latc+")")
  dlist = rotate(ugrid,vgrid,lonc,latc,0)
  ur = dlist[0]
  vr = dlist[1]
  printVarSummary(ur)
  printMinMax(ur,0)
  printVarSummary(vr)
  printMinMax(vr,0)
  ;; divide symmetric & asymmetric components
  uz = zonalAve(ur)
  vz = zonalAve(vr)
  us = ur
  us = conform(ur, uz, 0)
  us({:75},:) = 0.0 ;; only TC effects
  vs = vr
  vs = conform(vr, vz, 0)
  vs({:75},:) = 0.0 ;; only TC effects
  ua = ur
  ua = ur - us
  va = vr
  va = vr - vs
  us@long_name = "Symmetric tangential wind"
  vs@long_name = "Symmetric radial wind"
  ua@long_name = "Asymmetric tangential wind"
  va@long_name = "Asymmetric radial wind"
  printVarSummary(us)
  printMinMax(us,0)
  printVarSummary(vs)
  printMinMax(vs,0)
  printVarSummary(ua)
  printMinMax(ua,0)
  printVarSummary(va)
  printMinMax(va,0)
  ;; inverse rotate
  slist = rotate(us,vs,lonc,latc,1)
  alist = rotate(ua,va,lonc,latc,1)
  delete(us)
  delete(vs)
  delete(ua)
  delete(va)
  us = slist[0]
  vs = slist[1]
  ua = alist[0]
  va = alist[1]
  printVarSummary(us)
  printMinMax(us,0)
  printVarSummary(vs)
  printMinMax(vs,0)
  printVarSummary(ua)
  printMinMax(ua,0)
  printVarSummary(va)
  printMinMax(va,0)
;  exit

  wss = us
  wss = sqrt( us^2 + vs^2 )
  wss@long_name = "Wind speed"
  wss@units = "m/s"
  printVarSummary(wss)
  printMinMax(wss,0)
  wsa = ua
  wsa = sqrt( ua^2 + va^2 )
  wsa@long_name = "Wind speed"
  wsa@units = "m/s"
  printVarSummary(wsa)
  printMinMax(wsa,0)

  t0 = time(0)
  t1 = time(it)
  dt = t1
  dt = (t1 - t0)/60 ;;minutes=>hours
  ft = sprinti("%0.2d",dt)
  date = cd_string(t1,"%H00UTC, %D")

;
; These lines were commented out due to the fact that they make the plot
; unmanageably large.
; wtype@wkWidth = 2500
; wtype@wkHeight = 2500
  figdir = "Tl479L100"
  out = "stream"+clev+"_asym_"+init+"+"+ft+"h_"+prtb
  wks = gsn_open_wks(wtype,figdir+"/"+out)

  mpres                       = True
  mpres@gsnDraw = False
  mpres@gsnFrame = False
  mpres@gsnMaximize           = True

;---Zoom in on map
  mpres@mpMinLatF             = -10
  mpres@mpMaxLatF             =  50
  mpres@mpMinLonF             = 120
  mpres@mpMaxLonF             = 180

  mpres@mpFillOn              = False
  mpres@mpOutlineOn           = True
  mpres@mpDataBaseVersion     = "MediumRes"   ; slightly better outlines
  ;mpres@mpLandFillColor = "gray91"
;---Data is regional
  mpres@gsnAddCyclic          = False

  mpres@gsnLeftString = "@"+clev+" hPa"

;---Change the streamline levels
;  mpres@stLevelSelectionMode = "ExplicitLevels"
;  mpres@stLevels       = (/0./)
;  mpres@stLevelPalette = "temp_diff_1lev"
;  mpres@stLevelPalette        = "BlueRed" 
;  mpres@stLevelSelectionMode  = "ManualLevels"
;  mpres@stMinLevelValF        = -10.
;  mpres@stMaxLevelValF        = 10.
;  mpres@stLevelSpacingF       = 1.

  mpres@cnFillOn = True
  mpres@cnLinesOn = False
  mpres@cnLineLabelsOn = False
  mpres@cnLevelSelectionMode = "ManualLevels"
  if(plev1.lt.500) then
  mpres@cnMinLevelValF = 5.0
  mpres@cnMaxLevelValF = 40.0
  mpres@cnLevelSpacingF = 5.0
  else
  mpres@cnMinLevelValF = 3.
  mpres@cnMaxLevelValF = 27.0
  mpres@cnLevelSpacingF = 3.
  end if
  gsn_define_colormap(wks,"sunshine_9lev")
  mpres@gsnSpreadColors = True
;  mpres@lbLabelBarOn = False

  mpres@stLineThicknessF      = 2.0         ; default is 1.0
  mpres@stMinDistanceF = 0.1
  ;mpres@stMinLineSpacingF = 0.007
  mpres@stMinArrowSpacingF = 0.035
  mpres@stArrowLengthF = 0.005
  ;mpres@stLineColor = "orange"
  mpres@pmTickMarkDisplayMode = "Always"    ; more detailed tickmarks

;  plots = new(2,graphic)

;  mpres@gsnCenterString          = date + " sym"
;  plots(0) = gsn_csm_streamline_contour_map_ce(wks, us, vs, wss, mpres)

  mpres@gsnCenterString          = ""
  mpres@gsnRightString          = date + ": " + title ;+ " asym"
  if(plev0.ne.plev1) then
    mpres@gsnCenterStringParallelPosF = 0.6
  end if
  plot = gsn_csm_streamline_contour_map_ce(wks, ua, va, wsa, mpres)

  dum = new(2,graphic)
  mres = True
  mres@gsMarkerIndex  = 5
  mres@gsMarkerSizeF  = 12.0          ; marker size
  mres@gsMarkerThicknessF = 4.0
  mres@gsMarkerColor  = "red"      ; maker color
  ;do i=0,1
  ;  dum(i) = gsn_add_polymarker(wks,plots(i),lonc,latc,mres)
  ;end do
  markers = gsn_add_polymarker(wks,plot,lonc,latc,mres)

  draw(plot)
  frame(wks)
;  pres = True
;  pres@gsnPanelLabelBar = True
;  gsn_panel(wks,plots,(/1,2/),pres)
end
