function acos(x) {
  return atan2(sqrt(1 - x * x), x)
}

function dist(lon1, lat1, lon2, lat2) {
  a = 6.371e3
  pi = acos(-1)
  deg2rad = pi / 180
  dlon = (lon1 - lon2) * deg2rad
  lat1 *= deg2rad
  lat2 *= deg2rad
  return a * acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dlon))
}

BEGIN { ft = 0 }
{
  yyyymmddhh = sprintf("%0.4d%0.2d%0.2d%0.2d", $1, $2, $3, $4)
  if (FILENAME==ARGV[1]) {
    lon[yyyymmddhh] = $5
    lat[yyyymmddhh] = $6
  } else {
    if (lon[yyyymmddhh] && lat[yyyymmddhh]) {
      print ft, dist(lon[yyyymmddhh], lat[yyyymmddhh], $5, $6)
      ft += 6
    }
  }
}
