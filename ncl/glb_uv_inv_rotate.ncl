;*************************************************
; NCL tutorial script: vort_1.ncl
;*************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "../rotate/librotate.ncl"
load "write_nc.ncl"
load "calc_e.ncl"

procedure rotate_pole(u, v, lon, ur, vr, opt)
; u,v dimension is (time,lev,{89.5:90},lon) or (time,lev,{-90:-89.5},lon)
; North pole; x-axis=>90E, y-axis=>180  South pole; x-axis=>90E, y-axis=>0
local coslon, sinlon
begin
  coslon = cos(lon*deg2rad)
  sinlon = sin(lon*deg2rad)
  if(opt.eq.0) then
  ;North pole
    ur =  u * conform_dims(dimsizes(u), coslon, 3) + v * conform_dims(dimsizes(v), sinlon, 3)
    vr = -u * conform_dims(dimsizes(u), sinlon, 3) + v * conform_dims(dimsizes(v), coslon, 3)
  else if(opt.eq.1) then
  ;South pole
    ur = u * conform_dims(dimsizes(u), coslon, 3) - v * conform_dims(dimsizes(v), sinlon, 3)
    vr = u * conform_dims(dimsizes(u), sinlon, 3) + v * conform_dims(dimsizes(v), coslon, 3)
  end if
  end if
end

begin
;************************************************
; variable and file handling
;************************************************
  ;init = "2019100900"
  ;mem = "mean"
  yyyy = str_get_cols(init, 0, 3)
  mm   = str_get_cols(init, 4, 5)
  mmddhh = str_get_cols(init, 4, 9)
  fu  = addfile("../../netcdf/tigge/"+yyyy+"/jma/rotate/np_glb_"+init+"_"+mem+"_ug.nc","r")  ; open netcdf file
  fv  = addfile("../../netcdf/tigge/"+yyyy+"/jma/rotate/np_glb_"+init+"_"+mem+"_vg.nc","r")  ; open netcdf file
  fp  = addfile("../../netcdf/tigge/"+yyyy+"/jma/rotate/np_glb_"+init+"_"+mem+"_slp.nc","r")  ; open netcdf file
  fh  = addfile("../../netcdf/tigge/"+yyyy+"/jma/rotate/np_glb_"+init+"_"+mem+"_gh.nc","r")  ; open netcdf file
  trackfile = "../pytrack/jma/gtrack"+init+"_"+mem+".txt"
  nd = 21
  
  ur   = fu->var;(:,:,180:360,180:360)                       ; pull u off file
  vr   = fv->var;(:,:,180:360,180:360)                       ; pull v off file
  pr   = fp->var;(:,:,180:360,180:360)                       ; pull v off file
  hr   = fh->var;(:,:,180:360,180:360)                       ; pull v off file
  t0  = fu->time
  lev0 = fu->level
  lat0 = dble2flt(fu->lat)
  lon0 = dble2flt(fu->lon)
  uout = "../../netcdf/tigge/"+yyyy+"/jma/rotate/inv_r"+sprinti("%d", latr)+"_glb_"+init+"_"+mem+"_ug.nc"
  vout = "../../netcdf/tigge/"+yyyy+"/jma/rotate/inv_r"+sprinti("%d", latr)+"_glb_"+init+"_"+mem+"_vg.nc"
  pout = "../../netcdf/tigge/"+yyyy+"/jma/rotate/inv_r"+sprinti("%d", latr)+"_glb_"+init+"_"+mem+"_slp.nc"
  hout = "../../netcdf/tigge/"+yyyy+"/jma/rotate/inv_r"+sprinti("%d", latr)+"_glb_"+init+"_"+mem+"_gh.nc"
  
  if (isfilepresent(uout)) then
    system("rm " + uout)
  end if
  if (isfilepresent(vout)) then
    system("rm " + vout)
  end if
  if (isfilepresent(pout)) then
    system("rm " + pout)
  end if
  if (isfilepresent(hout)) then
    system("rm " + hout)
  end if
;************************************************
; prepare filtering function
;************************************************
  ;latr0 = 77.5
  latr0 = 90. - latr
  r0 = latr0
  deg2r(latr0, r0)
  print(r0)
  r = lat0({latr0:})
  lenr = dimsizes(r)
  deg2r(lat0({latr0:}), r)
  E = new(lenr(0), float)
  calc_e(r, r0, E)
  E = conform(E, 1.0, -1) - E
;************************************************
; remove symmetric wind around TC
;************************************************
  u = ur
  v = vr
  u_ = ur
  v_ = vr
  u_ = conform_dims(dimsizes(u), dim_avg_n(u, 3), (/0,1,2/))
  v_ = conform_dims(dimsizes(v), dim_avg_n(v, 3), (/0,1,2/))
  u(:, :, {latr0:}, :) = u(:, :, {latr0:}, :) - u_(:, :, {latr0:}, :) ;* conform_dims(dimsizes(u_(:, :, {latr0:}, :)), E, 2)
  v(:, :, {latr0:}, :) = v(:, :, {latr0:}, :) - v_(:, :, {latr0:}, :) ;* conform_dims(dimsizes(v_(:, :, {latr0:}, :)), E, 2)
  copy_VarCoords(ur, u)
  copy_VarCoords(vr, v)
  u@units = "m/s"
  v@units = "m/s"
  printMinMax(u, 0)
  printMinMax(v, 0)
  delete(ur)
  delete(vr)
;************************************************
; remove symmetric msl around TC
;************************************************
  p = pr
  p_ = pr
  pb = dim_avg_n(p(:,{latr0},:), 1)
  p_ = conform_dims(dimsizes(p), dim_avg_n(p, 2), (/0,1/)) - conform_dims(dimsizes(p), pb, 0)
  p(:, {latr0:}, :) = p(:, {latr0:}, :) - p_(:, {latr0:}, :) ;* conform_dims(dimsizes(p_(:, {latr0:}, :)), E, 1)
  copy_VarCoords(pr, p)
  p@units = "Pa"
  printMinMax(p, 0)
  delete(pr)
  ;exit
;************************************************
; remove symmetric gh around TC
;************************************************
  h = hr
  h_ = hr
  hb = dim_avg_n(h(:, :, {latr0}, :), 2)
  h_ = conform_dims(dimsizes(h), dim_avg_n(h, 3), (/0,1,2/)) - conform_dims(dimsizes(h), hb, (/0,1/))
  h(:, :, {latr0:}, :) = h(:, :, {latr0:}, :) - h_(:, :, {latr0:}, :) ;* conform_dims(dimsizes(p_(:, {latr0:}, :)), E, 1)
  copy_VarCoords(hr, h)
  h@units = "gpm"
  printMinMax(h, 0)
  delete(hr)
  ;exit
;************************************************
; rotate coordinate
;************************************************
  ndim = dimsizes(u)
  ntim = ndim(0)
  nlev = ndim(1)
  nlat = ndim(2)
  nlon = ndim(3)
  lonout = new(nlon * nlat, "float")
  latout = new(nlon * nlat, "float")
  track = asciiread(trackfile, (/nd, 7/), "float")
  ur = u
  vr = v
  pr = p
  hr = h
  do l = 0, nd-1
    dd = sprinti("%0.2d",toint(track(l, 2)))
    hh = sprinti("%0.2d",toint(track(l, 3)))
    print(dd+hh)
    lonc = track(l, 4) * deg2rad
    latc = track(l, 5) * deg2rad
    print("lonc = "+sprintf("%5.3f",lonc*rad2deg)+", latc = "+sprintf("%5.3f",latc*rad2deg))
    lonin = lon0 * deg2rad
    latin = lat0 * deg2rad
    inv_rotate_lonlat(lonc, latc, lonin, latin, lonout, latout)
    printMinMax(lonout, 0)
    printMinMax(latout, 0)
    uin = u(l,:,:,:)
    vin = v(l,:,:,:)
    pin = p(l,:,:)
    hin = h(l,:,:,:)
; rotate msl
    ptc = linint2_points(lon0, lat0, pin, True, lonout*rad2deg, latout*rad2deg, 0)
    if(any(ismissing(ptc)))
      print("Missing values in interpolation of msl")
      pmsg = linmsg(ptc(:), 0)
      ptc = (/pmsg/)
    end if
    pr(l,:,:) = (/reshape(ptc, (/nlat,nlon/))/)
; rotate gh
    htc = linint2_points(lon0, lat0, hin, True, lonout*rad2deg, latout*rad2deg, 0)
    if(any(ismissing(htc)))
      print("Missing values in interpolation of gh")
      do k = 0, nlev-1
        hmsg = linmsg(htc(k,:),0)
        htc(k,:) = (/hmsg/)
      end do
    end if
    do k = 0, nlev-1
      hr(l,k,:,:) = (/reshape(htc(k,:),(/nlat,nlon/))/)
    end do
; (u, v) => (xdot, ydot, zdot)
    xd = uin
    yd = uin
    zd = uin
    lon3d = conform(uin, lonin, 2)
    lat3d = conform(uin, latin, 1)
    uv2xyzd(uin, vin, lon3d, lat3d, xd, yd, zd) 
    if(any(ismissing(xd)).and.any(ismissing(yd)).and.any(ismissing(zd)))
      print("Missing values in xyzd")
      exit
    end if   
    xdnp = linint2_points(lon0, lat0, xd(:,:,:), True, lonout*rad2deg, latout*rad2deg, 0)
    ydnp = linint2_points(lon0, lat0, yd(:,:,:), True, lonout*rad2deg, latout*rad2deg, 0)
    zdnp = linint2_points(lon0, lat0, zd(:,:,:), True, lonout*rad2deg, latout*rad2deg, 0)
    if(any(ismissing(xdnp)).and.any(ismissing(ydnp)).and.any(ismissing(zdnp)))
      print("Missing values in interpolation")
      ;lmiss = new(nlev, "logical")
      do k = 0, nlev-1
        nsize = dimsizes(lonout)
        do ij = 0, nsize(0)-1
          if(ismissing(xdnp(k,ij)))
            ;lmiss(k) = True
            print(sprinti("%d",ij))
            print("level "+sprintf("%4.1f",lev0(k))+", lon = "+sprintf("%5.3f",lonout(ij)*rad2deg)+", lat = "+sprintf("%5.3f",latout(ij)*rad2deg))
          end if
        end do
      end do
      datamsg = linmsg(xdnp(:,:),0)
      xdnp = (/datamsg/)
      datamsg = linmsg(ydnp(:,:),0)
      ydnp = (/datamsg/)
      datamsg = linmsg(zdnp(:,:),0)
      zdnp = (/datamsg/)
      ;delete(datamsg)
      ;exit
    end if
; rotate coordinates to polar grid
    xdtc = xdnp
    ydtc = xdnp
    zdtc = xdnp
    lonc2d = conform(xdtc, lonc, -1)
    latc2d = conform(xdtc, latc, -1)
    np2tc(lonc2d, latc2d, xdnp, ydnp, zdnp, xdtc, ydtc, zdtc)
    if(any(ismissing(xdtc)).and.any(ismissing(ydtc)).and.any(ismissing(zdtc)))
      print("Missing values in rotation")
      exit
    end if
; (xdot, ydot, zdot) => (u, v)
    utc = xdtc
    vtc = xdtc
    lontc =conform(xdtc, ndtooned(conform_dims((/nlat, nlon/), lonin, 1)), 1)
    xyzd2uv(xdtc, ydtc, zdtc, lontc, utc, vtc)
    if(any(ismissing(utc)).and.any(ismissing(vtc)))
      print("Missing values in tc")
      exit
    end if
    ;printVarSummary(dataout)
    do k = 0, nlev-1
      ;ij = 0
      ;do j = 0, nlat-1
      ;  do i = 0, nlon-1
      ur(l,k,:,:) = (/reshape(utc(k,:),(/nlat,nlon/))/)
      vr(l,k,:,:) = (/reshape(vtc(k,:),(/nlat,nlon/))/)
      ;    ij = ij + 1
      ;  end do
      ;end do
    end do
  end do

;************************************************
; create netcdf4 file
;************************************************
  write_nc(uout, ur, t0, lev0, lat0, lon0)
  write_nc(vout, vr, t0, lev0, lat0, lon0)
  write_nc2d(pout, pr, t0, lat0, lon0)
  write_nc(hout, hr, t0, lev0, lat0, lon0)
  ;write_nc(ofile,vrt_r,t0,lev0,lat0,lon0)
  ;exit()

;************************************************
; create plot
;************************************************
  wks  = gsn_open_wks("png","uv_inv")              ; send graphics to PNG file

  res                   = True                   ; plot mods desired
  ;res@gsnPolar          = "SH"
  ;es@cnFillOn          = True                   ; turn on color
  ;res@cnFillPalette     = "BlWhRe"               ; set color map
  ;res@cnLinesOn         = False                  ; no contour lines

  ;res@gsnAddCyclic	      = False 
  res@mpLandFillColor      = "gray"            ; set land to be gray
  res@mpMinLonF            =  120.               ; select a subregion
  res@mpMaxLonF            =  160.
  res@mpMinLatF            =  10. 
  res@mpMaxLatF            =  30.

  ;res@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  ;res@cnMinLevelValF       = -20                   ; set min contour level
  ;res@cnMaxLevelValF       =  20                   ; set max contour level
  ;res@cnLevelSpacingF      =  2.                 ; set contour spacing

  res@gsnScalarContour      = True              ; vectors over scalars
  res@vcRefMagnitudeF       = 20.0              ; vec magnitude
  res@vcRefLengthF          = 0.05              ; size of ref vector
  res@vcLineArrowColor      = "black"           ; white vectors
  res@vcMinDistanceF        = 0.1              ; thins arrows near pole  
  res@vcGlyphStyle          = "CurlyVector"     ; turns on curly vectors
  res@vcLineArrowThicknessF = 2.
  res@vcMinMagnitudeF       = 5.0

  ;vrt = vrt*scale
  ;plot = gsn_csm_contour_map_polar(wks,vrt(3,:,:),res)
  plot = gsn_csm_vector_scalar_map(wks,ur(6,3,:,:),vr(6,3,:,:),hr(6,3,:,:),res)
  ;plot = gsn_csm_contour_map(wks,vrt(0,4,:,:),res) ; create plot  
;************************************************
end
