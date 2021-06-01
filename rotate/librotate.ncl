
pi = acos(-1.0)
tau = 2 * pi
deg2rad = pi / 180.0
rad2deg = 180.0 / pi

procedure lonlat2xyz(lon, lat, x, y, z)
begin
  x = cos(lon) * cos(lat)
  y = sin(lon) * cos(lat)
  z = sin(lat)
end

procedure xyz2lonlat(x, y, z, lon, lat)
begin
  lon = where(x .ne. 0.0, mod(atan2(y, x) + tau, tau), 0.0)
  lat= asin(z)
end

procedure uv2xyzd(u, v, lon, lat, xd, yd, zd)
local coslon, sinlon, sinlat
begin
  coslon = cos(lon)
  sinlon = sin(lon)
  sinlat = sin(lat)
  xd = -u * sinlon - v * coslon * sinlat
  yd =  u * coslon - v * sinlon * sinlat
  zd =               v * cos(lat)
end

procedure xyzd2uv(xd, yd, zd, lon, u, v)
local coslon, sinlon
begin
  coslon = cos(lon)
  sinlon = sin(lon)
  u = (/yd * coslon - xd * sinlon/)
  v = (/sign_f90(sqrt((xd * coslon + yd * sinlon)^2 + zd^2), zd)/)
end

procedure generate_points(nlon, nlat, dlat, lon, lat)
local dlon
begin
  dlon = 2 * pi / nlon
  lon = dlon * ispan(0, nlon - 1, 1)
  lat = 0.5 * pi - dlat * ispan(0, nlat - 1, 1)
end

procedure np2tc(lonc, latc, x, y, z, xx, yy, zz)
local coslonc, sinlonc, coslatc, sinlatc
begin
  coslonc = cos(lonc)
  sinlonc = sin(lonc)
  coslatc = cos(latc)
  sinlatc = sin(latc)
  xx = coslonc * sinlatc * x - sinlonc * y + coslonc * coslatc * z
  yy = sinlonc * sinlatc * x + coslonc * y + sinlonc * coslatc * z
  zz =          -coslatc * x               +           sinlatc * z
end

procedure tc2np(lonc, latc, x, y, z, xx, yy, zz)
local coslonc, sinlonc, coslatc, sinlatc
begin
  coslonc = cos(lonc)
  sinlonc = sin(lonc)
  coslatc = cos(latc)
  sinlatc = sin(latc)
  xx = coslonc * sinlatc * x + sinlonc * sinlatc * y - coslatc * z
  yy =          -sinlonc * x +           coslonc * y
  zz = coslonc * coslatc * x + sinlonc * coslatc * y + sinlatc * z
end

procedure rotate_lonlat(lonc, latc, lon, lat, lonout, latout)
local x, y, z, xx, yy, zz
begin
  nlon = dimsizes(lon)
  nlat = dimsizes(lat)
  ij = 0
  do j = 0, nlat(0) - 1
    do i = 0, nlon(0) - 1
      x = 0.0
      y = 0.0
      z = 0.0
      lonlat2xyz(lon(i), lat(j), x, y, z)
      xx = 0.0
      yy = 0.0
      zz = 0.0
      np2tc(lonc, latc, x, y, z, xx, yy, zz)
      xyz2lonlat(xx, yy, zz, lonout(ij), latout(ij))
      ij = ij + 1
    end do
  end do
end

procedure inv_rotate_lonlat(lonc, latc, lon, lat, lonout, latout)
local x, y, z, xx, yy, zz
begin
  nlon = dimsizes(lon)
  nlat = dimsizes(lat)
  ij = 0
  do j = 0, nlat(0) - 1
    do i = 0, nlon(0) - 1
      x = 0.0
      y = 0.0
      z = 0.0
      lonlat2xyz(lon(i), lat(j), x, y, z)
      xx = 0.0
      yy = 0.0
      zz = 0.0
      tc2np(lonc, latc, x, y, z, xx, yy, zz)
      xyz2lonlat(xx, yy, zz, lonout(ij), latout(ij))
      ij = ij + 1
    end do
  end do
end
