import sys
import netCDF4
import numpy as np
from datetime import datetime,timedelta
from pathlib import Path
import pandas as pd
import librotate

#Usage echo yyyymmddhh | python create_netcdf.py
#param = sys.stdin.readline().strip("\n").split(" ")
#yyyymmddhh = param[0]

datadir = Path('./')
init = sys.argv[1] #yyyymmddhh
yyyy = init[0:4]
mm   = init[4:6]
outdir = Path(f'/Volumes/dandelion/netcdf/era5/{yyyy}/{mm}')

nc_sfc = Path(f'np_sfc_{init}.nc')
nc_pl  = Path(f'np_pl_{init}.nc')
nc_uv  = Path(f'np_uv_{init}.nc')

outvar_sfc = ['msl','sp','t2m','d2m']
outvar_pl = ['t','phi','r','omg','pv','vo']
outvar_uv = ['u10', 'v10', 'u', 'v']
stname = {
'msl':'Pressure at sea level',
'sp':'Pressure at ground surface',
'tp':'Total precipitaion before 3 hours',
't2m':'Temparature at 2m-height',
'd2m':'Dew point temperature at 2m-height',
't':'Temparature',
'phi':'Geopotential',
'r':'Relative Humidity',
'omg':'Vertical P-Velocity',
'pv':'Potential Vorticity',
'vo':'Relative Vorticity',
'u10':'U component of wind at 10m-height',
'v10':'V component of wind at 10m-height', 
'u':'U component of wind',
'v':'V component of wind'}
units = {'msl':'Pa','sp':'Pa','tp':'m','t2m':'K',
'd2m':'K','t':'K','phi':'m^2/s^2','r':'%','omg':'Pa/s',
'pv':'K*m^2/kg/s','vo':'/s'}

outvar_dict = {}

start = datetime.strptime(init, '%Y%m%d%H')
dt = timedelta(hours=6)

outnc = netCDF4.Dataset(outdir/Path(f'np_{init}.nc'), 'w')
# create output information about dimensions and variables
in_sfc = netCDF4.Dataset(datadir/nc_sfc,'r')
nlat = in_sfc.variables["lat"][:].size
nlon = in_sfc.variables["lon"][:].size
print(f"nlat={nlat},nlon={nlon}")

time = outnc.createDimension('time',None)
level = outnc.createDimension('level',37)
lat = outnc.createDimension('lat',nlat)
lon = outnc.createDimension('lon',nlon)

times = outnc.createVariable("time",np.float64,("time",))
times.units = "seconds since 1970-01-01 00:00:00.0 0:00"
times.calendar = "standard"
times.axis = "T"
outvar_dict["time"] = times

levels = outnc.createVariable("level",np.float64,("level",))
levels.units = "hPa"
levels.axis = "Z"
outvar_dict["level"] = levels

latitudes = outnc.createVariable("lat",np.float64,("lat",))
latitudes.units = "degrees_north"
latitudes.axis = "Y"
outvar_dict["lat"] = latitudes

longitudes = outnc.createVariable("lon",np.float64,("lon",))
longitudes.units = "degrees_east"
longitudes.axis = "X"
outvar_dict["lon"] = longitudes

nvar_sfc = len(outvar_sfc)
nvar_pl = len(outvar_pl)
nvar_uv = len(outvar_uv)
for i in range(nvar_sfc):
    vkey = outvar_sfc[i]
    #if vkey == "T":
    #    vkey = "T2m"
    #elif vkey == "RH":
    #    vkey = "RH2m"
    var = outnc.createVariable(vkey,np.float32,("time","lat","lon"))
    var.standard_name = stname[vkey]
    var.units = units[vkey]
    outvar_dict[vkey] = var
for i in range(nvar_pl):
    vkey = outvar_pl[i]
    var = outnc.createVariable(vkey,np.float32,("time","level","lat","lon"))
    var.standard_name = stname[vkey]
    var.units = units[vkey]
    outvar_dict[vkey] = var
for i in range(nvar_uv):
    vkey = outvar_uv[i]
    if i < 2:
        var = outnc.createVariable(vkey,np.float32,("time","lat","lon"))
    else:
        var = outnc.createVariable(vkey,np.float32,("time","level","lat","lon"))
    var.standard_name = stname[vkey]
    var.units = "m/s"
    outvar_dict[vkey] = var
print(outvar_dict)

# surface
in_sfc = netCDF4.Dataset(datadir/nc_sfc,'r')
print(in_sfc.variables["time"].units)
#outvar_dict["level"][:] = in_scl.variables["level"][:]
outvar_dict["lat"][:] = in_sfc.variables["lat"][:]
outvar_dict["lon"][:] = in_sfc.variables["lon"][:]
for t in range(len(in_sfc.variables["time"][:])):
    date = netCDF4.num2date(in_sfc.variables["time"][t],in_sfc.variables["time"].units)
    t0 = netCDF4.date2num(date,outvar_dict["time"].units)
    print(date,t0)
    outvar_dict["time"][t] = t0
    for name in in_sfc.variables.keys():
        if(name == "time" or name == "lat" \
           or name == "lon" or name == "level"):
            continue
        elif(name in outvar_dict):
            print(name)
            outvar_dict[name][t,:] = in_sfc.variables[name][t,:]
            print(outvar_dict[name][:].shape)
        else:
            print("error")
# pressure levels
in_pl = netCDF4.Dataset(datadir/nc_pl,'r')
print(in_pl.variables["time"].units)
outvar_dict["level"][:] = in_pl.variables["level"][:]
for t in range(len(in_pl.variables["time"][:])):
    date = netCDF4.num2date(in_pl.variables["time"][t],in_pl.variables["time"].units)
    t0 = netCDF4.date2num(date,outvar_dict["time"].units)
    print(date,t0)
    outvar_dict["time"][t] = t0
    for name in in_pl.variables.keys():
        if(name == "time" or name == "lat" \
           or name == "lon" or name == "level"):
            continue
        if(name in outvar_dict):
            print(name)
            outvar_dict[name][t,:] = in_pl.variables[name][t,:]
            print(outvar_dict[name][:].shape)
        else:
            print("error")
# vector
in_uv = netCDF4.Dataset(datadir/nc_uv,'r')
for t in range(len(in_uv.variables["time"][:])):
    for name in in_uv.variables.keys():
        if(name == "time" or name == "lat" \
           or name == "lon" or name == "level"):
            continue
        if(name in outvar_dict):
            print(name)
            outvar_dict[name][t,:] = in_uv.variables[name][t,:]
            print(outvar_dict[name][:].shape)
        else:
            print("error")
print("time",outnc.variables["time"][:])
print("level",outnc.variables["level"][:])
print("lat",outnc.variables["lat"][:])
print("lon",outnc.variables["lon"][:])