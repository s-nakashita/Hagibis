load "librotate.ncl"
begin
  dlat = dlat * deg2rad
  lonin = new(nlon, "float")
  latin = new(nlat, "float")
  generate_points(nlon, nlat, dlat, lonin, latin)
  lonout = new(nlon * nlat , "float")
  latout = new(nlon * nlat , "float")

  system("rm -f "+outfile) 
  setfileoption("bin", "WriteByteOrder", "BigEndian")
  track = asciiread(trackfile, (/nstop/nfrpos + 1, 7/), "float")
  do i = 0, nstop, nfrpos
    r = i / nfrpos
    print("i="+i)
    lonc = track(r, 4) * deg2rad
    latc = track(r, 5) * deg2rad
    rotate_lonlat(lonc, latc, lonin, latin, lonout, latout)
    f = addfile(ncprefix+sprinti("%0.6d", i)+ncsuffix, "r")
    lon = doubletofloat(f->lon * deg2rad)
    lat = doubletofloat(f->lat * deg2rad)
    u = f->U(0,:,:,:)
    v = f->V(0,:,:,:)
    if(any(ismissing(u)).and.any(ismissing(v)))
      print("Missing values in input")
      exit
    end if
; (u, v) => (xdot, ydot, zdot)
    xd = u
    yd = u
    zd = u
    lon3d = conform(u, lon, 2)
    lat3d = conform(u, lat, 1)
    uv2xyzd(u, v, lon3d, lat3d, xd, yd, zd)
    if(any(ismissing(xd)).and.any(ismissing(yd)).and.any(ismissing(zd)))
      print("Missing values in xyzd")
      exit
    end if
; bilinear interpolation of winds
    xdtc = linint2_points(lon, lat(::-1), xd(:,::-1,:), True, lonout, latout, 0)
    ydtc = linint2_points(lon, lat(::-1), yd(:,::-1,:), True, lonout, latout, 0)
    zdtc = linint2_points(lon, lat(::-1), zd(:,::-1,:), True, lonout, latout, 0)
    if(any(ismissing(xdtc)).and.any(ismissing(ydtc)).and.any(ismissing(zdtc)))
      print("Missing values in interpolation")
      exit
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
    fbindirwrite(outfile, unp)
    fbindirwrite(outfile, vnp)
  end do
end
