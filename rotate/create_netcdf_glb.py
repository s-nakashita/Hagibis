import sys
import netCDF4
import numpy as np
from datetime import datetime,timedelta
from pathlib import Path
import logging
import pandas as pd
import librotate

logging.basicConfig(level=logging.DEBUG)
#Usage echo yyyymmddhh | python create_netcdf.py
param = sys.stdin.readline().strip("\n").split(" ")
yyyymmddhh = param[0]

datadir = Path('./')
outdir  = Path('./')
nc_scl  = Path('np_glb_sc_' + yyyymmddhh + '_mean.nc')
nc_vec  = Path('np_glb_ve_' + yyyymmddhh + '_mean.nc')
#nc_scl  = Path('np_sc_anl.nc')
#nc_vec  = Path('np_ve_anl.nc')

outvar = ['msl','t2m','d2m',\
          'u10','v10',\
          'u','v','t','q','gh']
#outvar = ['PRES_meansealevel','TMP_1D5maboveground','DPT_1D5maboveground',\
#          'UGRD_10maboveground','VGRD_10maboveground',\
#          'UGRD','VGRD','TMP','SPFH','HGT']
stname = ['Pressure','Temparature','Dew point of temparature',\
          'U component of wind','V component of wind',\
          'U component of wind','V component of wind',\
          'Temparature','Specific humidity','Geopotential height']
unit = ['Pa','K','K','m/s','m/s','m/s','m/s','K','kg/kg','m',]
outvar_dict = {}

start = datetime.strptime(yyyymmddhh, '%Y%m%d%H')
dt = timedelta(hours=6)

outnc = netCDF4.Dataset(outdir/Path('np_glb_'+yyyymmddhh+'_mean.nc'), 'w')
#outnc = netCDF4.Dataset(outdir/Path('np_anl.nc'), 'w')
in_scl = netCDF4.Dataset(datadir/nc_scl,'r')
nlat = in_scl.variables["latitude"][:].size
nlon = in_scl.variables["longitude"][:].size
logging.debug(f"nlat={nlat} nlon={nlon}")

time = outnc.createDimension('time',None)
level = outnc.createDimension('level',5)
lat = outnc.createDimension('lat',nlat)
lon = outnc.createDimension('lon',nlon)

times = outnc.createVariable("time",np.float64,("time",))
times.units = "hours since 1970-01-01 00:00:00.0 0:00"
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

for i in range(len(outvar)):
    vkey = outvar[i]
    if(i < 5):
        var = outnc.createVariable(vkey,np.float32,("time","lat","lon"))
    else:
        var = outnc.createVariable(vkey,np.float32,("time","level","lat","lon"))
    var.standard_name = stname[i]
    var.units = unit[i]
    outvar_dict[vkey] = var
logging.debug(outvar_dict)

in_scl = netCDF4.Dataset(datadir/nc_scl,'r')
logging.debug(in_scl.variables["time"].units)
logging.debug(in_scl.variables["level"][:])
outvar_dict["level"][:] = in_scl.variables["level"][:]
outvar_dict["lat"][:] = in_scl.variables["latitude"][:]
outvar_dict["lon"][:] = in_scl.variables["longitude"][:]
for t in range(len(in_scl.variables["time"][:])):
    date = netCDF4.num2date(in_scl.variables["time"][t],in_scl.variables["time"].units)
    t0 = netCDF4.date2num(date,outvar_dict["time"].units)
    #logging.debug(f"date {date}, t0 {t0}")
    outvar_dict["time"][t] = t0
    for name in in_scl.variables.keys():
        if(name == "time" or name == "latitude" \
           or name == "longitude" or name == "level"):
            continue
        if(name in outvar_dict):
            logging.info(name)
            outvar_dict[name][t,:] = in_scl.variables[name][t,:]
            logging.debug(outvar_dict[name][:].shape)
        else:
            logging.error("error")
            exit

in_vec = netCDF4.Dataset(datadir/nc_vec,'r')
logging.debug(in_vec.variables["level"][:])
for t in range(len(in_vec.variables["time"][:])):
    for name in in_vec.variables.keys():
        if(name == "time" or name == "latitude" \
           or name == "longitude" or name == "level"):
            continue
        if(name in outvar_dict):
            logging.info(name)
            outvar_dict[name][t,:] = in_vec.variables[name][t,:]
            logging.debug(outvar_dict[name][:].shape)
        else:
            logging.error("error")
            exit

logging.info(outnc.variables["time"][:])
logging.info(outnc.variables["level"][:])
logging.info(outnc.variables["lat"][:])
logging.info(outnc.variables["lon"][:])
