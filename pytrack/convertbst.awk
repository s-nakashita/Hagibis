{
  y = 20 substr($1, 1, 2)
  m = substr($1, 3, 2) * 1
  d = substr($1, 5, 2) * 1
  h = substr($1, 7, 2) * 1
  print y, m, d, h, $5 * 0.1, $4*0.1, $6 * 100
}
