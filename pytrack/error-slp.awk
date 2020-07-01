function rmse(slp1, slp2) {
  slp1 *= 0.01
  slp2 *= 0.01
  return sqrt((slp1 - slp2) * (slp1 - slp2))
}

BEGIN { ft = 0 }
{
  yyyymmddhh = sprintf("%0.4d%0.2d%0.2d%0.2d", $1, $2, $3, $4)
  if (FILENAME==ARGV[1]) {
    slp[yyyymmddhh] = $7
  } else {
    if (slp[yyyymmddhh]) {
      print ft, rmse(slp[yyyymmddhh], $7)
      ft += 6
    }
  }
}
