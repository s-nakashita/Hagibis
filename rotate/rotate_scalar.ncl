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
    lon = f->lon
    lat = f->lat
    datain = f->$varname$
    if (dimsizes(dimsizes(datain)) .eq. 4) then
      dataout = linint2_points(lon, lat(::-1), datain(0,:,::-1,:), True, lonout * rad2deg, latout * rad2deg, 0)
    else
      dataout = linint2_points(lon, lat(::-1), datain(0,::-1,:), True, lonout * rad2deg, latout * rad2deg, 0)
    end if
    if(any(ismissing(dataout)))
      print("Missing values in interpolation")
    end if
    fbindirwrite(outfile, dataout)
  end do
end
