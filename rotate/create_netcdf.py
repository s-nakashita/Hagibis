import sys
import netCDF4
import numpy as np
from datetime import datetime,timedelta
from pathlib import Path
import pandas as pd
import librotate

#Usage echo yyyymmddhh | python create_netcdf.py
param = sys.stdin.readline().strip("\n").split(" ")
yyyymmddhh = param[0]

datadir = Path('./')
outdir  = Path('./')
nc_scl  = Path('np_sc_' + yyyymmddhh + '_mean.nc')
nc_vec  = Path('np_ve_' + yyyymmddhh + '_mean.nc')
#nc_scl  = Path('np_sc_anl.nc')
#nc_vec  = Path('np_ve_anl.nc')

outvar = ['PRES_meansealevel','TMP_2maboveground','DPT_2maboveground',\
          'UGRD_10maboveground','VGRD_10maboveground',\
          'UGRD','VGRD','TMP','SPFH','HGT']
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

outnc = netCDF4.Dataset(outdir/Path('np_'+yyyymmddhh+'_mean.nc'), 'w')
#outnc = netCDF4.Dataset(outdir/Path('np_anl.nc'), 'w')
in_scl = netCDF4.Dataset(datadir/nc_scl,'r')
nlat = in_scl.variables["latitude"][:].size
nlon = in_scl.variables["longitude"][:].size
print(nlat,nlon)

time = outnc.createDimension('time',None)
level = outnc.createDimension('level',5)
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

for i in range(len(outvar)):
    vkey = outvar[i]
    if(i < 5):
        var = outnc.createVariable(vkey,np.float32,("time","lat","lon"))
    else:
        var = outnc.createVariable(vkey,np.float32,("time","level","lat","lon"))
    var.standard_name = stname[i]
    var.units = unit[i]
    outvar_dict[vkey] = var
print(outvar_dict)

in_scl = netCDF4.Dataset(datadir/nc_scl,'r')
print(in_scl.variables["time"].units)
outvar_dict["level"][:] = in_scl.variables["level"][:]
outvar_dict["lat"][:] = in_scl.variables["latitude"][:]
outvar_dict["lon"][:] = in_scl.variables["longitude"][:]
for t in range(len(in_scl.variables["time"][:])):
    date = netCDF4.num2date(in_scl.variables["time"][t],in_scl.variables["time"].units)
    t0 = netCDF4.date2num(date,outvar_dict["time"].units)
    print(date,t0)
    outvar_dict["time"][t] = t0
    for name in in_scl.variables.keys():
        if(name == "time" or name == "latitude" \
           or name == "longitude" or name == "level"):
            continue
        if(name in outvar_dict):
            print(name)
            outvar_dict[name][t,:] = in_scl.variables[name][t,:]
            print(outvar_dict[name][:].shape)
        else:
            print(int(name[-5:-2]))
            lev=1
            if(int(name[-5:-2]) < 925):
                lev+=1
            if(int(name[-5:-2]) < 850):
                lev+=1
            if(int(name[-5:-2]) < 700):
                lev+=1
            if(int(name[-5:-2]) < 500):
                lev+=1
            if(int(name[-5:-2]) < 300):
                lev+=1
            if(int(name[-5:-2]) < 250):
                lev+=1
            if(int(name[-5:-2]) == 0):
                lev=0
                        
            if(name[:3]=="TMP"):
                print(name[:3])
                outvar_dict["TMP"][t,lev,:] = in_scl.variables[name][t,:]
                print(outvar_dict[name[:3]][t])
            elif(name[:4]=="SPFH"):
                print(name[:4])
                outvar_dict["SPFH"][t,lev,:] = in_scl.variables[name][t,:]
                print(outvar_dict[name[:4]][t])
            elif(name[:3]=="HGT"):
                print(name[:3])
                outvar_dict["HGT"][t,lev,:] = in_scl.variables[name][t,:]
                print(outvar_dict[name[:3]][t])
            else:
                print("error")

in_vec = netCDF4.Dataset(datadir/nc_vec,'r')
for t in range(len(in_vec.variables["time"][:])):
    for name in in_vec.variables.keys():
        if(name == "time" or name == "latitude" \
           or name == "longitude" or name == "level"):
            continue
        if(name in outvar_dict):
            print(name)
            outvar_dict[name][t,:] = in_vec.variables[name][t,:]
            print(outvar_dict[name][:].shape)
        else:
            print(int(name[-5:-2]))
            lev=1
            if(int(name[-5:-2]) < 925):
                lev+=1
            if(int(name[-5:-2]) < 850):
                lev+=1
            if(int(name[-5:-2]) < 700):
                lev+=1
            if(int(name[-5:-2]) < 500):
                lev+=1
            if(int(name[-5:-2]) < 300):
                lev+=1
            if(int(name[-5:-2]) < 250):
                lev+=1
            if(int(name[-5:-2]) == 0):
                lev=0
                        
            if(name[:4]=="UGRD"):
                print(name[:4])
                outvar_dict["UGRD"][t,lev,:] = in_vec.variables[name][t,:]
                print(outvar_dict[name[:4]][t])
            elif(name[:4]=="VGRD"):
                print(name[:4])
                outvar_dict["VGRD"][t,lev,:] = in_vec.variables[name][t,:]
                print(outvar_dict[name[:4]][t])
            else:
                print("error")

print(outnc.variables["time"][:],outnc.variables["lat"][:],outnc.variables["lon"][:])
