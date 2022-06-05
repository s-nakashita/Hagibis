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
exp = sys.argv[1]
init = sys.argv[2] #yyyymmddhh
if exp == "cntl":
    outdir = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est')
elif exp == "prtb1":
    outdir = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p')
elif exp == "prtb2":
    outdir = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p2')
elif exp == "prtbn":
    outdir = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+pn')
elif exp == "prtbf":
    outdir = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+pf')
else:
    logger.error('don\'t exist such a experiment')
    exit()
nc_sfc = Path(f'np_fcst_surf_{init}.nc')
nc_pl  = Path(f'np_fcst_p_{init}.nc')
nc_uv  = Path(f'np_fcst_uv_{init}.nc')

outvar_sfc = ['PHI','PSEA','P','RAIN','T2m','RH2m','CLA','CLL','CLM','CLH']
outvar_pl = ['T','Z','RH','OMG','CWC','CVR','VOR']
outvar_uv = ['U10m', 'V10m', 'U', 'V']
#outvar = ['PRES_meansealevel','TMP_1D5maboveground','DPT_1D5maboveground',\
#          'UGRD_10maboveground','VGRD_10maboveground',\
#          'UGRD','VGRD','TMP','SPFH','HGT']
stname = ['Geopotential','Pressure at sea level','Pressure','Rain','Temparature','Relative Humidity',\
    'Cloud Fraction', 'Cloud Fraction (Low)', 'Cloud Fraction (Middle)', 'Cloud Fraction (High)',\
    'Temparature','Geopotential height','Relative Humidity','Vertical Velocity(Pressure)','Cloud Water Content','Cloud Cover','Relative Vorticity',\
    'U component of wind', 'V component of wind', 'U component of wind', 'V component of wind']
unit = ['kgm^2/s^2','Pa','Pa','mm','K','percent',\
    'm^2/m^2','m^2/m^2','m^2/m^2','m^2/m^2',\
    'K','m','percent','Pa/s','kg/kg','m^2/m^2','/s',\
    'm/s','m/s','m/s','m/s']
outvar_dict = {}

start = datetime.strptime(init, '%Y%m%d%H')
dt = timedelta(hours=6)

outnc = netCDF4.Dataset(outdir/Path(f'np_fcst_{init}.nc'), 'w')
# create output information about dimensions and variables
in_sfc = netCDF4.Dataset(datadir/nc_sfc,'r')
nlat = in_sfc.variables["lat"][:].size
nlon = in_sfc.variables["lon"][:].size
print(f"nlat={nlat},nlon={nlon}")

time = outnc.createDimension('time',None)
level = outnc.createDimension('level',22) #global:22, asia:21
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
    var.standard_name = stname[i]
    var.units = unit[i]
    outvar_dict[vkey] = var
for i in range(nvar_pl):
    vkey = outvar_pl[i]
    var = outnc.createVariable(vkey,np.float32,("time","level","lat","lon"))
    var.standard_name = stname[i+nvar_sfc]
    var.units = unit[i+nvar_sfc]
    outvar_dict[vkey] = var
for i in range(nvar_uv):
    vkey = outvar_uv[i]
    if i < 2:
        var = outnc.createVariable(vkey,np.float32,("time","lat","lon"))
    else:
        var = outnc.createVariable(vkey,np.float32,("time","level","lat","lon"))
    var.standard_name = stname[i+nvar_sfc+nvar_pl]
    var.units = unit[i+nvar_sfc+nvar_pl]
    outvar_dict[vkey] = var
print(outvar_dict)

# surface
in_sfc = netCDF4.Dataset(datadir/nc_sfc,'r')
print(in_sfc.variables["time"].units)
#outvar_dict["level"][:] = in_scl.variables["level"][:]
outvar_dict["lat"][:] = in_sfc.variables["lat"][:]
outvar_dict["lon"][:] = in_sfc.variables["lon"][:]
t_sfc = []
for t in range(len(in_sfc.variables["time"][:])):
    date = netCDF4.num2date(in_sfc.variables["time"][t],in_sfc.variables["time"].units)
    t0 = netCDF4.date2num(date,outvar_dict["time"].units)
    print(date,t0)
    outvar_dict["time"][t] = t0
    t_sfc.append(t0)
    for name in in_sfc.variables.keys():
        if(name == "time" or name == "lat" \
           or name == "lon" or name == "level"):
            continue
        if name == "T":
            print(name)
            outvar_dict["T2m"][t,:] = in_sfc.variables[name][t,:]
            print(outvar_dict["T2m"][:].shape)
        elif name == "RH":
            print(name)
            outvar_dict["RH2m"][t,:] = in_sfc.variables[name][t,:]
            print(outvar_dict["RH2m"][:].shape)
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
    #outvar_dict["time"][t] = t0
    index = t_sfc.index(t0)
    for name in in_pl.variables.keys():
        if(name == "time" or name == "lat" \
           or name == "lon" or name == "level"):
            continue
        if(name in outvar_dict):
            print(name)
            outvar_dict[name][index,:] = in_pl.variables[name][t,:]
            print(outvar_dict[name][:].shape)
        else:
            print("error")
# vector
in_uv = netCDF4.Dataset(datadir/nc_uv,'r')
for t in range(len(in_uv.variables["time"][:])):
    date = netCDF4.num2date(in_uv.variables["time"][t],in_uv.variables["time"].units)
    t0 = netCDF4.date2num(date,outvar_dict["time"].units)
    print(date,t0)
    #outvar_dict["time"][t] = t0
    index = t_sfc.index(t0)
    for name in in_uv.variables.keys():
        if(name == "time" or name == "lat" \
           or name == "lon" or name == "level"):
            continue
        if(name in outvar_dict):
            print(name)
            outvar_dict[name][index,:] = in_uv.variables[name][t,:]
            print(outvar_dict[name][:].shape)
        else:
            print("error")
print("time",outnc.variables["time"][:])
print("level",outnc.variables["level"][:])
print("lat",outnc.variables["lat"][:])
print("lon",outnc.variables["lon"][:])