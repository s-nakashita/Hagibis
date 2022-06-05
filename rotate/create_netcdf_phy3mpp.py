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
if exp == "clim":
    outdir = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl959L100')
elif exp == "est":
    outdir = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est')
elif exp == "mgdsst":
    outdir = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst')
else:
    logger.error('don\'t exist such a experiment')
    exit()
nc_sfc = Path(f'np_fcst_phy2m_avr_{init}.nc')
nc_pl  = Path(f'np_fcst_phy3mpp_{init}.nc')
#nc_uv  = Path(f'np_fcst_uv_asia_{init}.nc')

outvar_sfc = ['PRCh','PRLh','FLSH','FLLH',
'RSDT','RSDB','RSUT','RSUB','RLDB','RLUT','RLUB']
outvar_pl = ['HRRS','HRRL','HRCV','HRLC','HRVD','HRAD','HR',
'QRCV','QRLC','QRVD','QRAD']
#outvar_uv = ['U10m', 'V10m', 'U', 'V']
#outvar = ['PRES_meansealevel','TMP_1D5maboveground','DPT_1D5maboveground',\
#          'UGRD_10maboveground','VGRD_10maboveground',\
#          'UGRD','VGRD','TMP','SPFH','HGT']
stname = {
'PRCh':'Precipitation (Convective)',
'PRLh':'Precipitation (Stratiform)',
'FLSH':'Specific Heat Flux at Ground Surface',
'FLLH':'Latent Heat Flux ad Ground Surface',
'RSDT':'Radiation Flux (Short,Down,Top)',
'RSDB':'Radiation Flux (Short,Down,Bottom)',
'RSUT':'Radiation Flux (Short,Up,Top)',
'RSUB':'Radiation Flux (Short,Up,Bottom)',
'RLDB':'Radiation Flux (Long,Down,Bottom)',
'RLUT':'Radiation Flux (Long,Up,Top)',
'RLUB':'Radiation Flux (Long,Up,Bottom)',
'HRRS':'Heating Rate (Radiative, short)',
'HRRL':'Heating Rate (Radiative, long)',
'HRCV':'Heating Rate (Convective)',
'HRLC':'Heating Rate (Stratiform)',
'HRVD':'Heating Rate (Turbulent)',
'HRAD':'Heating Rate (Advective)',
'HR':'Heating Rate (Total)',
'QRCV':'Moistening Rate (Convective)',
'QRLC':'Moistening Rate (Stratiform)',
'QRVD':'Moistening Rate (Turbulent)',
'QRAD':'Moistening Rate (Advective)'}
unit = {'PRCh':'kg/m^2/h','PRLh':'kg/m^2/h',
'FLSH':'W/m^2','FLLH':'W/m^2','RSDT':'W/m^2',
'RSDB':'W/m^2','RSUT':'W/m^2','RSUB':'W/m^2',
'RLDB':'W/m^2','RLUT':'W/m^2','RLUB':'W/m^2',
'HRRS':'K/h','HRRL':'K/h','HRCV':'K/h',
'HRLC':'K/h','HRVD':'K/h','HRAD':'K/h','HR':'K/h',
'QRCV':'kg/kg/h','QRLC':'kg/kg/h','QRVD':'kg/kg/h','QRAD':'kg/kg/h'}
outvar_dict = {}

start = datetime.strptime(init, '%Y%m%d%H')
dt = timedelta(hours=6)

outnc = netCDF4.Dataset(outdir/Path(f'np_fcst_phy_{init}.nc'), 'w')
# create output information about dimensions and variables
in_pl = netCDF4.Dataset(datadir/nc_pl,'r')
nlat = in_pl.variables["lat"][:].size
nlon = in_pl.variables["lon"][:].size
print(f"nlat={nlat},nlon={nlon}")

time = outnc.createDimension('time',None)
level = outnc.createDimension('level',21)
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
#nvar_uv = len(outvar_uv)
for i in range(nvar_sfc):
    vkey = outvar_sfc[i]
    if vkey == "PRCd":
        vkey = "PRCh"
    elif vkey == "PRLd":
        vkey = "PRLh"
    var = outnc.createVariable(vkey,np.float32,("time","lat","lon"))
    var.standard_name = stname[vkey]
    var.units = unit[vkey]
    outvar_dict[vkey] = var
for i in range(nvar_pl):
    vkey = outvar_pl[i]
    var = outnc.createVariable(vkey,np.float32,("time","level","lat","lon"))
    var.standard_name = stname[vkey]
    var.units = unit[vkey]
    outvar_dict[vkey] = var
#for i in range(nvar_uv):
#    vkey = outvar_uv[i]
#    if i < 2:
#        var = outnc.createVariable(vkey,np.float32,("time","lat","lon"))
#    else:
#        var = outnc.createVariable(vkey,np.float32,("time","level","lat","lon"))
#    var.standard_name = stname[i+nvar_sfc+nvar_pl]
#    var.units = unit[i+nvar_sfc+nvar_pl]
#    outvar_dict[vkey] = var
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
        if name == "PRCd":
            print(name)
            outvar_dict["PRCh"][t,:] = in_sfc.variables[name][t,:]
            print(outvar_dict["PRCh"][:].shape)
        elif name == "PRLd":
            print(name)
            outvar_dict["PRLh"][t,:] = in_sfc.variables[name][t,:]
            print(outvar_dict["PRLh"][:].shape)
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
outvar_dict["lat"][:] = in_pl.variables["lat"][:]
outvar_dict["lon"][:] = in_pl.variables["lon"][:]
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
## vector
#in_uv = netCDF4.Dataset(datadir/nc_uv,'r')
#for t in range(len(in_uv.variables["time"][:])):
#    for name in in_uv.variables.keys():
#        if(name == "time" or name == "lat" \
#           or name == "lon" or name == "level"):
#            continue
#        if(name in outvar_dict):
#            print(name)
#            outvar_dict[name][t,:] = in_uv.variables[name][t,:]
#            print(outvar_dict[name][:].shape)
#        else:
#            print("error")
print("time",outnc.variables["time"][:])
print("level",outnc.variables["level"][:])
print("lat",outnc.variables["lat"][:])
print("lon",outnc.variables["lon"][:])