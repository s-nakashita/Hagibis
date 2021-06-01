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
  ;mem="mean"
  yyyy = str_get_cols(init, 0, 3)
  mm   = str_get_cols(init, 4, 5)
  mmddhh = str_get_cols(init, 4, 9)
  ;in  = addfile("../../netcdf/tigge/"+yyyy+"/jma/glb_"+mmddhh+"_mean.nc","r")  ; open netcdf file
  in  = addfile("../../netcdf/tigge/"+yyyy+"/jma/glb_"+init+"_"+mem+".nc","r")  ; open netcdf file
  u   = in->u;(:,:,180:360,180:360)                       ; pull u off file
  v   = in->v;(:,:,180:360,180:360)                       ; pull v off file
  msl = in->msl
  gh  = in->gh
  t0  = in->time
  lev0 = in->level
  lat0 = in->latitude
  lon0 = in->longitude
  trackfile = "../pytrack/jma/gtrack"+init+"_"+mem+".txt"
  nd = 21
  
  uout = "../../netcdf/tigge/"+yyyy+"/jma/rotate/np_glb_"+init+"_"+mem+"_ug.nc"
  vout = "../../netcdf/tigge/"+yyyy+"/jma/rotate/np_glb_"+init+"_"+mem+"_vg.nc"
  pout = "../../netcdf/tigge/"+yyyy+"/jma/rotate/np_glb_"+init+"_"+mem+"_slp.nc"
  hout = "../../netcdf/tigge/"+yyyy+"/jma/rotate/np_glb_"+init+"_"+mem+"_gh.nc"
  
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
; rotate pole vector
;************************************************
  print(dim_avg_n(u(0,0,{89.5},:),0))
  print(dim_avg_n(v(0,0,{89.5},:),0))
  print(dim_avg_n(u(0,0,{90},:),0))
  print(dim_avg_n(v(0,0,{90},:),0))
  unp = u(:,:,{89.5:90},:)
  vnp = v(:,:,{89.5:90},:)
  ur = new(dimsizes(unp),typeof(unp))
  vr = new(dimsizes(vnp),typeof(vnp))
  rotate_pole(unp, vnp, lon0, ur, vr, 0)
  print(dim_avg_n(ur(0,0,0,:),0))
  u(:,:,{89.5:90},:) = (/ur/)
  v(:,:,{89.5:90},:) = (/vr/)
  usp = u(:,:,{-90:-89.5},:)
  vsp = v(:,:,{-90:-89.5},:)
  rotate_pole(usp, vsp, lon0, ur, vr, 1)
  print(dim_avg_n(ur(0,0,0,:),0))
  u(:,:,{-90:-89.5},:) = (/ur/)
  v(:,:,{-90:-89.5},:) = (/vr/)
  ;print(u(0,0,{-90},360))
  ;print(v(0,0,{-90},360))
  delete(unp)
  delete(vnp)
  delete(ur)
  delete(vr)
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
  pr = msl
  hr = gh
  do l = 0, nd-1
    dd = sprinti("%0.2d",toint(track(l, 2)))
    hh = sprinti("%0.2d",toint(track(l, 3)))
    print(dd+hh)
    lonc = track(l, 4) * deg2rad
    latc = track(l, 5) * deg2rad
    print("lonc = "+sprintf("%5.3f",lonc*rad2deg)+", latc = "+sprintf("%5.3f",latc*rad2deg))
    lonin = lon0 * deg2rad
    latin = lat0 * deg2rad
    rotate_lonlat(lonc, latc, lonin, latin, lonout, latout)
    printMinMax(lonout, 0)
    printMinMax(latout, 0)
    uin = u(l,:,:,:)
    vin = v(l,:,:,:)
    pin = msl(l,:,:)
    hin = gh(l,:,:,:)
; rotate msl
    pnp = linint2_points(lon0, lat0, pin, True, lonout*rad2deg, latout*rad2deg, 0)
    if(any(ismissing(pnp)))
      print("Missing values in interpolation of msl")
      pmsg = linmsg(pnp(:),0)
      pnp = (/pmsg/)
    end if
    pr(l,:,:) = (/reshape(pnp,(/nlat,nlon/))/)
; rotate msl
    hnp = linint2_points(lon0, lat0, hin, True, lonout*rad2deg, latout*rad2deg, 0)
    if(any(ismissing(hnp)))
      print("Missing values in interpolation of gh")
      do k = 0, nlev-1
        hmsg = linmsg(hnp(k,:),0)
        hnp(k,:) = (/hmsg/)
      end do
    end if
    do k = 0, nlev-1
      hr(l,k,:,:) = (/reshape(hnp(k,:),(/nlat,nlon/))/)
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
    xdtc = linint2_points(lon0, lat0, xd(:,:,:), True, lonout*rad2deg, latout*rad2deg, 0)
    ydtc = linint2_points(lon0, lat0, yd(:,:,:), True, lonout*rad2deg, latout*rad2deg, 0)
    zdtc = linint2_points(lon0, lat0, zd(:,:,:), True, lonout*rad2deg, latout*rad2deg, 0)
    if(any(ismissing(xdtc)).and.any(ismissing(ydtc)).and.any(ismissing(zdtc)))
      print("Missing values in interpolation")
      ;lmiss = new(nlev, "logical")
      do k = 0, nlev-1
        nsize = dimsizes(lonout)
        do ij = 0, nsize(0)-1
          if(ismissing(xdtc(k,ij)))
            ;lmiss(k) = True
            print(sprinti("%d",ij))
            print("level "+sprintf("%4.1f",lev0(k))+", lon = "+sprintf("%5.3f",lonout(ij)*rad2deg)+", lat = "+sprintf("%5.3f",latout(ij)*rad2deg))
          end if
        end do
      end do
      datamsg = linmsg(xdtc(:,:),0)
      xdtc = (/datamsg/)
      datamsg = linmsg(ydtc(:,:),0)
      ydtc = (/datamsg/)
      datamsg = linmsg(zdtc(:,:),0)
      zdtc = (/datamsg/)
      ;delete(datamsg)
      ;exit
    end if
; rotate coordinates to polar grid
    xdnp = xdtc
    ydnp = xdtc
    zdnp = xdtc
    lonc2d = conform(xdtc, lonc, -1)
    latc2d = conform(xdtc, latc, -1)
    tc2np(lonc2d, latc2d, xdtc, ydtc, zdtc, xdnp, ydnp, zdnp)
    if(any(ismissing(xdnp)).and.any(ismissing(ydnp)).and.any(ismissing(zdnp)))
      print("Missing values in rotation")
      exit
    end if
; (xdot, ydot, zdot) => (u, v)
    unp = xdnp
    vnp = xdnp
    lonnp =conform(xdnp, ndtooned(conform_dims((/nlat, nlon/), lonin, 1)), 1)
    xyzd2uv(xdnp, ydnp, zdnp, lonnp, unp, vnp)
    if(any(ismissing(unp)).and.any(ismissing(vnp)))
      print("Missing values in np")
      exit
    end if
    ;printVarSummary(dataout)
    do k = 0, nlev-1
      ;ij = 0
      ;do j = 0, nlat-1
      ;  do i = 0, nlon-1
      ur(l,k,:,:) = (/reshape(unp(k,:),(/nlat,nlon/))/)
      vr(l,k,:,:) = (/reshape(vnp(k,:),(/nlat,nlon/))/)
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
  wks  = gsn_open_wks("png","uv")              ; send graphics to PNG file

  res                   = True                   ; plot mods desired
  res@gsnPolar          = "NH"
  res@mpMinLatF         = 80
  ;res@mpMaxLatF         = -80
  res@cnFillOn          = False                  ; turn on color
  res@cnFillPalette     = "BlWhRe"               ; set color map
  res@cnLinesOn         = True                  ; no contour lines

  ;res@gsnAddCyclic	      = False 
  ;res@mpLandFillColor      = "gray"            ; set land to be gray
  ;res@mpMinLonF            =  120.               ; select a subregion
  ;res@mpMaxLonF            =  160.
  ;res@mpMinLatF            =  10. 
  ;res@mpMaxLatF            =  60.

  ;res@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  ;res@cnMinLevelValF       = -20                   ; set min contour level
  ;res@cnMaxLevelValF       =  20                   ; set max contour level
  ;res@cnLevelSpacingF      =  2.                 ; set contour spacing

  res@gsnScalarContour      = True              ; vectors over scalars
  res@vcRefMagnitudeF       = 20.0              ; vec magnitude
  res@vcRefLengthF          = 0.05              ; size of ref vector
  res@vcLineArrowColor      = "black"           ; white vectors
  res@vcMinDistanceF        = 0.1              ; thins arrows near pole  
  res@vcGlyphStyle          = "LineArrow"     ; turns on curly vectors
  res@vcLineArrowThicknessF = 2.
  res@vcMinMagnitudeF       = 5.0

  ;vrt = vrt*scale
  ;plot = gsn_csm_contour_map_polar(wks,vrt(3,:,:),res)
  plot = gsn_csm_vector_scalar_map_polar(wks,ur(6,3,:,:),vr(6,3,:,:),hr(6,3,:,:),res)
  ;plot = gsn_csm_contour_map(wks,vrt(0,4,:,:),res) ; create plot  
;************************************************
end
