# rotate program for ERA5
import sys
import os
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import logging
import xarray as xr
import pandas as pd
import librotate

logging.basicConfig(level=logging.DEBUG)

nlon       = 1440 # every 0.25 degree
nlat       = 113  # every 0.25 degree
latmax     = 28.0 #degree
dlat = latmax/(nlat-1) #degree

init = sys.argv[1] #yyyymmddhh
yyyy = init[0:4]
mm   = init[4:6]
datadir = Path(f'/Volumes/dandelion/netcdf/era5/{yyyy}/{mm}')
trackf = Path(f'../pytrack/track_era5.txt')

outsfcnc   = f'np_sfc_{init}.nc'
outplnc   = f'np_pl_{init}.nc'
outdir = Path('./')
logging.debug(os.path.isfile(trackf))

var_sfc = ['msl','sp','t2m','d2m']
var_pl = ['t','z','r','w','pv','vo']

lonin,latin = librotate.generate_points(nlon,nlat,dlat)
logging.info(f"input coordinate lon:{lonin} \n lat:{latin}")
#print(lonin)
#print(latin)

sfc_np=[]
pl_np=[]
nt = 0
with trackf.open() as track:
    for l in track:
        l_split = l.split()
        logging.info(l_split)
        year  = l_split[0]
        month = l_split[1]
        day   = l_split[2]
        hour  = l_split[3]
        if day == "13" and hour == "12":
            break
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
        newlon = xr.DataArray(lonout,dims='np_lonlat')
        newlat = xr.DataArray(latout,dims='np_lonlat')
        ## record
        lonlatout = np.concatenate([lonout.reshape(-1,1), latout.reshape(-1,1)], 1)
        np.savetxt(
            f"outline_era5/lonlatout_{int(year):04d}{int(month):02d}{int(day):02d}{int(hour):02d}.txt", 
            lonlatout[-nlon:,:])
        
        insfcnc = f'{int(day):02d}/slev_{int(day):02d}.nc'
        inplnc  = f'{int(day):02d}/plev_{int(day):02d}.nc'
        
        # surface
        datain = xr.open_dataset(datadir / insfcnc).sel(time=date)
        logging.debug(datain)

        lon = datain["longitude"].values
        lat = datain["latitude"].values
        #lev = data["level"].values
        #nlev = len(lev)
        #print(lon,lat)
        ## add cyclic
        #databc = data.sel(lon=0.0)
        #databc["lon"] = 360.0
        #logging.debug(databc)
        #datain = xr.concat([data, databc], dim="lon")
        #logging.debug(datain)
        sfcout=[]
        for vname in var_sfc:
            din = datain[vname]
            logging.debug(din)
            ##missing values need to be searched
            df = din.values
            if(np.isnan(df).sum() != 0):
                logging.warning("missing value exists in input data.")
                continue
            ## (for RAIN) accumulate value => 3h accumulated value
            if vname == 'tp' and nt > 1:
                y_pre = int(year)
                m_pre = int(month)
                d_pre = int(day)
                h_pre = int(hour) - 3
                if h_pre < 0:
                    d_pre = d_pre - 1
                    h_pre = h_pre + 24
                date_str  = f'{y_pre:04d}/{m_pre:02d}/{d_pre:02d} {h_pre:02d}:00'
                logging.debug("3 hours before : "+date_str)
                date_pre = datetime.strptime(date_str,'%Y/%m/%d %H:%M')
                data_pre = xr.open_dataset(datadir / insfcnc).sel(time=date_pre)
                din = din - data_pre[vname]
            inattrs = din.attrs
            #print(inattrs)
            data_interp = din.interp(longitude=newlon, latitude=newlat)
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
                                    ('lat',latin),('lon',lonin)],\
                                   attrs=inattrs,name=vname)
            logging.debug(dataout)
            sfcout.append(dataout)
        # pressure levels
        datain = xr.open_dataset(datadir / inplnc).sel(time=date)
        logging.debug(datain)

        lon = datain["longitude"].values
        lat = datain["latitude"].values
        lev = datain["level"].values
        nlev = len(lev)
        #print(lon,lat)
        plout = []
        for vname in var_pl:
            #for lev in ['200','250','300','500','850']:
            #    vname_pl = vname.replace('000',lev)
            din = datain[vname]
            logging.debug(din)
            ##missing value
            df = din.values
            if(np.isnan(df).sum() != 0):
                logging.warning("missing value exists in input data.")
                continue
            inattrs = din.attrs
            #print(inattrs)
            data_interp = \
                din.interp(longitude=newlon, latitude=newlat)
            ##missing value
            #df = data_interp.to_pandas()
            value = data_interp.values
            if(np.isnan(value).sum() != 0):
                logging.warning("#{} missing value exists in interpolation data.".format\
                (np.isnan(value).sum()))
                msk = np.isnan(value[0])
                mislon = newlon[msk]
                mislat = newlat[msk]
                for i in range(len(mislon)):
                    logging.warning("missing at lon{} lat{}".format(mislon[i].values,mislat[i].values))
                continue
            if vname == 'z':
                vname = 'phi'
            if vname == 'w':
                vname = 'omg'
            dataout = xr.DataArray(value.reshape(1,nlev,nlat,nlon),\
                                   [('time',pd.date_range(date,periods=1)),\
                                    ('level',lev),\
                                 ('lat',latin),('lon',lonin)],\
                                       attrs=inattrs,name=vname)
            logging.debug(dataout)
            plout.append(dataout)

        sfc_np.append(xr.merge(sfcout))
        pl_np.append(xr.merge(plout))
        nt += 1
        #date_pre = date
data_sfc_np = xr.merge(sfc_np)
logging.info(data_sfc_np)
data_sfc_np.to_netcdf(outdir / outsfcnc, 'w')
data_pl_np = xr.merge(pl_np)
logging.info(data_pl_np)
data_pl_np.to_netcdf(outdir / outplnc, 'w')