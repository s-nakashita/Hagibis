import sys
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import xarray as xr
import pandas as pd
import librotate

#Usage echo yyyymmddhh datadir trackf nlon nlat latmax | python rotate_scalar.py
param = sys.stdin.readline().strip("\n").split(" ")
yyyymmddhh = param[0]
ddirname   = param[1]
trackname  = param[2]
nlon       = int(param[3])
nlat       = int(param[4])
latmax     = float(param[5]) #degree
dlat = latmax/nlat #degree

datadir = Path(ddirname)
outdir = Path('./')
mmddhh = yyyymmddhh[4:]
innc = mmddhh+'_mean.nc'
outnc = 'np_sc_'+yyyymmddhh+'_mean.nc'
trackf = Path(trackname)

var_sfc = ['PRES_meansealevel','TMP_2maboveground','DPT_2maboveground']
#var_sfc = ['PRES_meansealevel','TMP_1D5maboveground','DPT_1D5maboveground']
#var_pl = ['TMP_000mb','SPFH_000mb','HGT_000mb']
var_pl = ['TMP','SPFH','HGT']

lonin,latin = librotate.generate_points(nlon,nlat,dlat)
print(lonin)
print(latin)

da_np=[]
with trackf.open() as track:
    for l in track:
        l_split = l.split()
        print(l_split)
        year  = l_split[0]
        month = l_split[1]
        day   = l_split[2]
        hour  = l_split[3]
        date_str  = year+'/'+month+'/'+day+' '+hour+':00'
        print(date_str)
        date = datetime.strptime(date_str,'%Y/%m/%d %H:%M')
        print(date)
        lonc  = float(l_split[4])
        latc  = float(l_split[5])
        print(lonc,latc)
    
        lonout,latout = librotate.rotate_lonlat(lonc,latc,lonin,latin)
        print(lonout.max(),lonout.min(),latout.max(),latout.min())
    
        data = xr.open_dataset(datadir / innc).sel(time=date)
        print(data)

        lon = data["lon"].values
        lat = data["lat"].values
        lev = data["level"].values
        nlev = len(lev)
        #print(lon,lat)
        newlon = xr.DataArray(lonout,dims='np_lonlat')
        newlat = xr.DataArray(latout,dims='np_lonlat')
        vname = var_sfc[0]
        daout=[]
        for vname in var_sfc:
            datain = data[vname]
            print(datain)
            ##missing values need to be searched
            df = datain.to_pandas()
            if(df.isnull().values.sum() != 0):
                print("missing value exists in input data.")
                continue
            inattrs = datain.attrs
            #print(inattrs)
            data_interp = datain.interp(lon=newlon, lat=newlat)
            #print(data_interp)
            ##missing value
            df = data_interp.to_pandas()
            if(df.isnull().values.sum() != 0):
                print("missing value exists in interpolation data.")
                continue
            value = data_interp.values
            #print(value)
            dataout = xr.DataArray(value.reshape(1,nlat,nlon),\
                                   [('time',pd.date_range(date,periods=1)),\
                                    ('latitude',latin),('longitude',lonin)],\
                                   attrs=inattrs,name=vname)
            print(dataout)
            daout.append(dataout)
            
        for vname in var_pl:
            #for lev in ['200','250','300','500','850']:
            #    vname_pl = vname.replace('000',lev)
            datain = data[vname]
            print(datain)
            ##missing value
            df = datain[0,:].to_pandas()
            if(df.isnull().values.sum() != 0):
                print("missing value exists in input data.")
                continue
            inattrs = datain.attrs
            #print(inattrs)
            data_interp = \
                datain.interp(lon=newlon, lat=newlat)
            ##missing value
            df = data_interp.to_pandas()
            if(df.isnull().values.sum() != 0):
                print("missing value exists in interpolation data.")
                continue
            value = data_interp.values
            dataout = xr.DataArray(value.reshape(1,nlev,nlat,nlon),\
                                   [('time',pd.date_range(date,periods=1)),\
                                    ('level',lev),\
                                 ('latitude',latin),('longitude',lonin)],\
                                       attrs=inattrs,name=vname)
            print(dataout)
            daout.append(dataout)

        da_np.append(xr.merge(daout))
data_np = xr.merge(da_np)
print(data_np)

data_np.to_netcdf(outdir / outnc, 'w')
