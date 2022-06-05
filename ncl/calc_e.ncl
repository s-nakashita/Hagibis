re = 6.381e3 ;km
pi = acos(-1.0)
tau = 2 * pi
deg2rad = pi / 180.0
rad2deg = 180.0 / pi
procedure deg2r(deg, r)
begin
  r = re * (conform(deg, 90.0, -1) - deg) * deg2rad
end
procedure r2deg(r, deg)
begin
  deg = conform(deg, 90.0, -1) - r / re * rad2deg 
end
procedure calc_e(r, r0, E)
begin
  l = 0.2*r0
  E = (exp(-((conform(r, r0, -1) - r)^2)/l/l) - conform(r, exp(-(r0^2)/l/l), -1)) / (1.0 - exp(-(r0^2)/l/l))
end