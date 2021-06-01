import sys
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import logging
import xarray as xr
import pandas as pd
import librotate

logging.basicConfig(level=logging.DEBUG)

#Usage echo yyyymmddhh datadir trackf nlon nlat latmax | python rotate_scalar.py
param = sys.stdin.readline().strip("\n").split(" ")
ddirname   = param[0]
trackname  = param[1]
nlon       = int(param[2])
nlat       = int(param[3])
latmax     = float(param[4]) #degree
dlat = latmax/(nlat-1) #degree

datadir = Path(ddirname)
outdir = Path('./')
innc = 'init.nc'
outnc = 'np_sc_init.nc'
trackf = Path(trackname)
logging.debug(trackf)

var_sfc = ['PRMSL_meansealevel','PRES_surface','TMP_2maboveground','RH_2maboveground',\
            'LCDC_surface','MCDC_surface','HCDC_surface','TCDC_surface']
#var_sfc = ['PRES_meansealevel','TMP_1D5maboveground','DPT_1D5maboveground']
#var_pl = ['TMP_000mb','SPFH_000mb','HGT_000mb']
var_pl = ['TMP','RH','HGT','VVEL']

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

        lon = data["lon"].values
        lat = data["lat"].values
        lev = data["level"].values
        nlev = len(lev)
        #print(lon,lat)
        newlon = xr.DataArray(lonout,dims='np_lonlat')
        newlat = xr.DataArray(latout,dims='np_lonlat')
        
        # add cyclic
        databc = data.sel(lon=0.0)
        databc["lon"] = 360.0
        logging.debug(databc)
        datain = xr.concat([data, databc], dim="lon")
        logging.debug(datain)
        daout=[]
        for vname in var_sfc:
            din = datain[vname]
            logging.debug(din)
            ##missing values need to be searched
            df = din.values
            if(np.isnan(df).sum() != 0):
                logging.warning("missing value exists in input data.")
                continue
            inattrs = din.attrs
            #print(inattrs)
            data_interp = din.interp(lon=newlon, lat=newlat)
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
            logging.debug(dataout)
            daout.append(dataout)
            
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
                din.interp(lon=newlon, lat=newlat)
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
