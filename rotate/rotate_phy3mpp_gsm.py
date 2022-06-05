# rotate program for GSM Tl959 experiments
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

exp = sys.argv[1]
init = sys.argv[2] #yyyymmddhh
if exp == "clim":
    datadir = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl959L100')
    trackf = Path(f'../pytrack/track{init}_gsm_tl959.txt')
elif exp == "est":
    datadir = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est')
    trackf = Path(f'../pytrack/track{init}_gsm_tl959_est.txt')
elif exp == "mgdsst":
    datadir = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst')
    trackf = Path(f'../pytrack/track{init}_gsm_tl959_mgdsst.txt')
else:
    logger.error('don\'t exist such a experiment')
    exit()
outdir = Path('./')
insfcnc = f'fcst_phy2m_avr_{init}.nc'
inplnc = f'fcst_phy3mpp_{init}.nc'
outsfcnc = 'np_'+insfcnc
outplnc = 'np_'+inplnc
logging.debug(os.path.isfile(trackf))

var_sfc = ['aPRCd','aPRLd','aFLSH','aFLLH',
'aRSDT','aRSDB','aRSUT','aRSUB','aRLDB','aRLUT','aRLUB']
var_pl = ['aHRRS','aHRRL','aHRCV','aHRLC','aHRVD','aHRAD','aHR',
'aQRCV','aQRLC','aQRVD','aQRAD']
units = {'aPRCd':'kg/m^2/h','aPRLd':'kg/m^2/h','aFLSH':'W/m^2','aFLLH':'W/m^2',
'aRSDT':'W/m^2','aRSDB':'W/m^2','aRSUT':'W/m^2','aRSUB':'W/m^2',
'aRLDB':'W/m^2','aRLUT':'W/m^2','aRLUB':'W/m^2',
'aHRRS':'K/h','aHRRL':'K/h','aHRCV':'K/h',
'aHRLC':'K/h','aHRVD':'K/h','aHRAD':'K/h','aHR':'K/h',
'aQRCV':'kg/kg/h','aQRLC':'kg/kg/h','aQRVD':'kg/kg/h','aQRAD':'kg/kg/h'}
lonin,latin = librotate.generate_points(nlon,nlat,dlat)
logging.info(f"input coordinate lon:{lonin} \n lat:{latin}")
#print(lonin)
#print(latin)
perday = 24.0 * 60.0 * 60.0 # secs
perhour = 60.0 * 60.0 # secs
factor = perhour / perday
sfc_np=[]
pl_np=[]
nt = 0
with trackf.open() as track:
    for l in track:
        if nt == 0:
            nt+=1
            continue
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
        newlon = xr.DataArray(lonout,dims='np_lonlat')
        newlat = xr.DataArray(latout,dims='np_lonlat')
        ## record
        #lonlatout = np.concatenate([lonout.reshape(-1,1), latout.reshape(-1,1)], 1)
        #np.savetxt(
        #    f"outline_{exp}/{init}/lonlatout_{exp}_{int(year):04d}{int(month):02d}{int(day):02d}{int(hour):02d}.txt", 
        #    lonlatout[-nlon:,:])
        
        # surface
        try:
            datain = xr.open_dataset(datadir / insfcnc).sel(time=date)
        except KeyError:
            logging.warning(f"{date} isn't in {insfcnc}")
            continue
        logging.debug(datain)

        lon = datain["lon"].values
        lat = datain["lat"].values
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
            ##(for PRCd and PRLd) accumulate value per day => 6h accumulated value per hour
            if vname == 'aPRCd' or vname == 'aPRLd':
                din = din * factor
            inattrs = din.attrs
            inattrs["units"] = units[vname]
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
            if vname == 'T' or vname == 'RH':
                vname = vname + "2m"
            dataout = xr.DataArray(value.reshape(1,nlat,nlon),\
                                   [('time',pd.date_range(date,periods=1)),\
                                    ('lat',latin),('lon',lonin)],\
                                   attrs=inattrs,name=vname[1:])
            logging.debug(dataout)
            sfcout.append(dataout)
        sfc_np.append(xr.merge(sfcout))
        # pressure levels
        try:
            datain = xr.open_dataset(datadir / inplnc).sel(time=date)
        except KeyError:
            logging.warning(f"{date} isn't in {inplnc}")
            continue
        logging.debug(datain)

        lon = datain["lon"].values
        lat = datain["lat"].values
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
            ## accumulate value per day => 6h accumulated value per hour
            #if vname[0] == 'a' and nt > 1:
            #    y_pre = int(year)
            #    m_pre = int(month)
            #    d_pre = int(day)
            #    h_pre = int(hour) - 6
            #    if h_pre < 0:
            #        d_pre = d_pre - 1
            #        h_pre = h_pre + 24
            #    date_str  = f'{y_pre:04d}/{m_pre:02d}/{d_pre:02d} {h_pre:02d}:00'
            #    logging.debug("6 hours before : "+date_str)
            #    date_pre = datetime.strptime(date_str,'%Y/%m/%d %H:%M')
            #    data_pre = xr.open_dataset(datadir / inplnc).sel(time=date_pre)
            #    din = (din - data_pre[vname])*factor
            din = din*factor
            inattrs = din.attrs
            inattrs["units"] = units[vname]
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
                                 ('lat',latin),('lon',lonin)],\
                                       attrs=inattrs,name=vname[1:])
            logging.debug(dataout)
            plout.append(dataout)

        pl_np.append(xr.merge(plout))
        nt += 1
        #date_pre = date
data_sfc_np = xr.merge(sfc_np)
logging.info(data_sfc_np)
data_sfc_np.to_netcdf(outdir / outsfcnc, 'w')
data_pl_np = xr.merge(pl_np)
logging.info(data_pl_np)
data_pl_np.to_netcdf(outdir / outplnc, 'w')