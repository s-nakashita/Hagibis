import sys
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import logging
import xarray as xr
import pandas as pd
import librotate

logging.basicConfig(level=logging.INFO)

#Usage echo yyyymmddhh datadir trackf nlon nlat latmax | python rotate_scalar.py
param = sys.stdin.readline().strip("\n").split(" ")
yyyymmddhh = param[0]
ddirname   = param[1]
trackname  = param[2]
nlon       = int(param[3])
nlat       = int(param[4])
latmax     = float(param[5]) #degree
dlat = latmax/(nlat-1) #degree

datadir = Path(ddirname)
outdir = Path('./')
mmddhh = yyyymmddhh[4:]
innc = 'glb_'+yyyymmddhh+'_mean.nc'
outnc = 'np_glb_sc_'+yyyymmddhh+'_mean.nc'
#innc = 'anl.nc'
#outnc = 'np_sc_anl.nc'
trackf = Path(trackname)
logging.debug(trackf)

var_sfc = ['msl','t2m','d2m']
var_pl = ['t','q','gh']

lonin,latin = librotate.generate_points(nlon,nlat,dlat)
logging.info(f"input coordinate lon:{lonin} \n lat:{latin}")
#print(lonin)
#print(latin)

da_np=[]
with trackf.open() as track:
    for l in track:
        l_split = l.split()
        logging.info(l_split)
        year  = l_split[0]
        month = l_split[1]
        day   = l_split[2]
        hour  = l_split[3]
        date_str  = year+'/'+month+'/'+day+' '+hour+':00'
        logging.debug(date_str)
        date = datetime.strptime(date_str,'%Y/%m/%d %H:%M')
        logging.info(date)
        lonc  = float(l_split[4])
        latc  = float(l_split[5])
        logging.info(f"TC center ({lonc},{latc})")
    
        lonout,latout = librotate.rotate_lonlat(lonc,latc,lonin,latin)
        logging.debug(\
            f"output coordinate lon:max{lonout.max():.3f}, min{lonout.min():.3f} \n \
                 lat:max{latout.max():.3f},min{latout.min():.3f}")
    
        data = xr.open_dataset(datadir / innc).sel(time=date)
        logging.debug(data)

        lon = data["longitude"].values
        lat = data["latitude"].values
        lev = data["level"].values
        nlev = len(lev)
        #print(lon,lat)
        newlon = xr.DataArray(lonout,dims='np_lonlat')
        newlat = xr.DataArray(latout,dims='np_lonlat')
        vname = var_sfc[0]
        daout=[]
        for vname in var_sfc:
            datain = data[vname]
            logging.debug(datain)
            ##missing values need to be searched
            df = datain.to_pandas()
            if(df.isnull().values.sum() != 0):
                logging.warning("missing value exists in input data.")
                continue
            inattrs = datain.attrs
            #print(inattrs)
            data_interp = datain.interp(longitude=newlon, latitude=newlat)
            logging.debug(data_interp)
            ##missing value
            #df = data_interp.to_pandas()
            value = data_interp.values
            if(np.isnan(value).sum() != 0):
                logging.warning("#{} missing value exists in interpolation data.".format\
                (np.isnan(value).sum()))
                msk = np.isnan(value)
                mislon = newlon[msk]
                mislat = newlat[msk]
                for i in range(len(mislon)):
                    logging.warning("missing at lon{} lat{}".format(mislon[i].values,mislat[i].values))
                continue
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
            logging.debug(datain)
            ##missing value
            df = datain[0,:].to_pandas()
            if(df.isnull().values.sum() != 0):
                logging.warning("missing value exists in input data.")
                continue
            inattrs = datain.attrs
            #print(inattrs)
            data_interp = \
                datain.interp(longitude=newlon, latitude=newlat)
            ##missing value
            #df = data_interp.to_pandas()
            value = data_interp.values
            if(np.isnan(value).sum() != 0):
                logging.warning("missing value exists in interpolation data.")
                continue
            dataout = xr.DataArray(value.reshape(1,nlev,nlat,nlon),\
                                   [('time',pd.date_range(date,periods=1)),\
                                    ('level',lev),\
                                 ('latitude',latin),('longitude',lonin)],\
                                       attrs=inattrs,name=vname)
            logging.debug(dataout)
            daout.append(dataout)

        da_np.append(xr.merge(daout))
data_np = xr.merge(da_np)
logging.info(data_np)

data_np.to_netcdf(outdir / outnc, 'w')
