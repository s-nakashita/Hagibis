;
; Write netcdf4 file script
;
procedure write_nc(ofile,var,time,lev,lat,lon)
; input: ofile(outfile name), var, time, lev, lat, lon(dims)
local out, ndim, ntim, nlev, nlat, nlon, dim_names, dim_sizes, dim_unlimited, tatts, levatts, latatts, lonatts, vatts
begin
out = addfile(ofile,"c")
ndim = dimsizes(var)
ntim = ndim(0)
nlev = ndim(1)
nlat = ndim(2)
nlon = ndim(3)
dim_names = (/"time","level","lat","lon"/)
dim_sizes = (/ntim,nlev,nlat,nlon/)
dim_unlimited = (/True,False,False,False/) ; time dimension is set unlimited
filedimdef(out,dim_names,dim_sizes,dim_unlimited)

print("Defining file dimensions")
filevardef(out,"time","double","time")
filevardef(out,"level","double","level")
filevardef(out,"lat","double","lat")
filevardef(out,"lon","double","lon")
  
print("Defining file dimension attributes")
tatts = True
tatts@standard_name = "time";t0@standard_name
tatts@units = time@units
filevarattdef(out,"time",tatts)

levatts = True
levatts@long_name = "pressure levels";lev0@long_name
levatts@units = lev@units
filevarattdef(out,"level",levatts)

latatts = True
latatts@long_name = "latitude";lat0@long_name
latatts@units = lat@units
filevarattdef(out,"lat",latatts)

lonatts = True
lonatts@long_name = "longitude";lon0@long_name
lonatts@units = lon@units
filevarattdef(out,"lon",lonatts)
out->level = (/lev/)
out->lat   = (/lat/)
out->lon   = (/lon/)

print("Defining file variable"+tochar(0))
filevardef(out,"var","float",(/"time","level","lat","lon"/))
vatts = True
vatts@long_name = var@long_name
vatts@nits      = var@units
vatts@_Fillvalue = 1e32
filevarattdef(out,"var",vatts)
do n = 0, ntim-1
    out->time(n) = (/time(n)/)
    out->var(n,:,:,:) = (/var(n,:,:,:)/)
end do

end

procedure write_nc2d(ofile,var,time,lat,lon)
; input: ofile(outfile name), var, time, lat, lon(dims)
local out, ndim, ntim, nlat, nlon, dim_names, dim_sizes, dim_unlimited, tatts, latatts, lonatts, vatts
begin
out = addfile(ofile,"c")
ndim = dimsizes(var)
ntim = ndim(0)
nlat = ndim(1)
nlon = ndim(2)
dim_names = (/"time","lat","lon"/)
dim_sizes = (/ntim,nlat,nlon/)
dim_unlimited = (/True,False,False/) ; time dimension is set unlimited
filedimdef(out,dim_names,dim_sizes,dim_unlimited)

print("Defining file dimensions")
filevardef(out,"time","double","time")
filevardef(out,"lat","double","lat")
filevardef(out,"lon","double","lon")
  
print("Defining file dimension attributes")
tatts = True
tatts@standard_name = "time";t0@standard_name
tatts@units = time@units
filevarattdef(out,"time",tatts)

latatts = True
latatts@long_name = "latitude";lat0@long_name
latatts@units = lat@units
filevarattdef(out,"lat",latatts)

lonatts = True
lonatts@long_name = "longitude";lon0@long_name
lonatts@units = lon@units
filevarattdef(out,"lon",lonatts)
out->lat   = (/lat/)
out->lon   = (/lon/)

print("Defining file variable"+tochar(0))
filevardef(out,"var","float",(/"time","lat","lon"/))
vatts = True
vatts@long_name = var@long_name
vatts@nits      = var@units
vatts@_Fillvalue = 1e32
filevarattdef(out,"var",vatts)
do n = 0, ntim-1
    out->time(n) = (/time(n)/)
    out->var(n,:,:) = (/var(n,:,:)/)
end do

end