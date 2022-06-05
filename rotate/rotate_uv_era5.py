import sys
import os
import numpy as np
from pathlib import Path
from datetime import datetime,timedelta
import logging
import xarray as xr
import pandas as pd
import librotate

logging.basicConfig(level=logging.DEBUG)

nlon       = 1440
nlat       = 113
latmax     = 28.0 #degree
dlat = latmax/(nlat-1) #degree

init = sys.argv[1] #yyyymmddhh
yyyy = init[0:4]
mm   = init[4:6]
datadir = Path(f'/Volumes/dandelion/netcdf/era5/{yyyy}/{mm}')
trackf = Path(f'../pytrack/track_era5.txt')

outnc   = f'np_uv_{init}.nc'
outdir = Path('./')
logging.debug(os.path.isfile(trackf))

var_sfc = ['u10','v10']
var_pl  = ['u','v']

lonin,latin = librotate.generate_points(nlon,nlat,dlat)
logging.info(f"input coordinate lon:{lonin} \n lat:{latin}")

da_np = []
with trackf.open() as track:
    for l in track:
        l_split = l.split()
        logging.info(l_split)
        year     = l_split[0]
        month    = l_split[1]
        day      = l_split[2]
        hour     = l_split[3]
        date_str = year + '/' + month + '/' + day + ' ' + hour + ':00'
        date     = datetime.strptime(date_str, '%Y/%m/%d %H:%M')
        logging.info(date)
        lonc     = float(l_split[4])
        latc     = float(l_split[5])
        logging.info(f"TC center ({lonc},{latc})")
        
        lonout,latout = librotate.rotate_lonlat(lonc,latc,lonin,latin)
        logging.debug(\
            f"output coordinate lon:max{lonout.max():.3f}, min{lonout.min():.3f} \n lat:max{latout.max():.3f},min{latout.min():.3f}")
        newlon = xr.DataArray(lonout,dims='np_lonlat')
        newlat = xr.DataArray(latout,dims='np_lonlat')

        insfcnc = f'{int(day):02d}/slev_{int(day):02d}.nc'
        inplnc  = f'{int(day):02d}/plev_{int(day):02d}.nc'
        
        daout = []
        #sfc
        datain = xr.open_dataset(datadir / insfcnc).sel(time=date)
        logging.debug(datain)

        lon = datain["longitude"].values
        lat = datain["latitude"].values
        #lev = data["level"].values
        #nlev = len(lev)
        #logging.debug(f"level: {lev}")
        u = datain[var_sfc[0]].values
        v = datain[var_sfc[1]].values 
        attrs_u = datain[var_sfc[0]].attrs
        attrs_v = datain[var_sfc[1]].attrs
        #missing values need to be searched
        #dfu = data[var_sfc[0]].to_pandas()
        #dfv = data[var_sfc[1]].to_pandas()
        if(np.isnan(u).sum() != 0 or np.isnan(v).sum() != 0):
            logging.warning("missing value exists in input")
            continue
        lon2d = lon.reshape(1,-1)
        lat2d = lat.reshape(-1,1)
        xd,yd,zd = librotate.uv2xyzd(u,v,lon2d,lat2d)
        #debug
        ud, vd = librotate.xyzd2uv(xd, yd, zd, lon2d, lat2d)
        logging.debug(f"uv2xyzd u {np.max(np.abs(u-ud))}")
        logging.debug(f"uv2xyzd v {np.max(np.abs(v-vd))}")
        da_xd = xr.DataArray(xd.reshape(1,len(lat),len(lon)),\
                         [('time',pd.date_range(date,periods=1)),('latitude',lat),('longitude',lon)],\
                         name='xd')
        da_yd = xr.DataArray(yd.reshape(1,len(lat),len(lon)),\
                         [('time',pd.date_range(date,periods=1)),('latitude',lat),('longitude',lon)],\
                         name='yd')
        da_zd = xr.DataArray(zd.reshape(1,len(lat),len(lon)),\
                         [('time',pd.date_range(date,periods=1)),('latitude',lat),('longitude',lon)],\
                         name='zd')
        
        xd_interp = da_xd.interp(longitude=newlon, latitude=newlat)
        yd_interp = da_yd.interp(longitude=newlon, latitude=newlat)
        zd_interp = da_zd.interp(longitude=newlon, latitude=newlat)
        #missing value
        dfx = xd_interp.values
        dfy = yd_interp.values
        dfz = zd_interp.values
        if(np.isnan(dfx).sum() != 0 or np.isnan(dfy).sum() != 0\
            or np.isnan(dfz).sum() != 0):
            logging.warning("missing value exists in interpolation")
            msk = np.isnan(dfx[0,:])
            mislon = newlon.values[msk]
            mislat = newlat.values[msk]
            for i in range(len(mislon)):
                logging.warning("missing at lon{} lat{}".format(mislon[i],mislat[i]))
            continue
        xdtc = xd_interp.values.reshape(nlat,nlon)
        ydtc = yd_interp.values.reshape(nlat,nlon)
        zdtc = zd_interp.values.reshape(nlat,nlon)
        xdnp,ydnp,zdnp = librotate.tc2np(lonc,latc,xdtc,ydtc,zdtc)
        #debug
        xb,yb,zb = librotate.np2tc(lonc,latc,xdnp,ydnp,zdnp)
        logging.debug(f"tc2np x {np.max(np.abs(xdtc-xb))}")
        logging.debug(f"tc2np y {np.max(np.abs(ydtc-yb))}")
        logging.debug(f"tc2np z {np.max(np.abs(zdtc-zb))}")
        #missing value
        if(np.isnan(xdnp).sum() != 0 or np.isnan(ydnp).sum() != 0 or \
            np.isnan(zdnp).sum() != 0):
            logging.warning("missing value exists in tc2np")
            continue
        lonnp = lonin.reshape(1,-1)
        latnp = latin.reshape(-1,1)
        #unp,vnp = librotate.xyzd2uv(xdnp,ydnp,zdnp,lonnp)
        unp,vnp = librotate.xyzd2uv(xdnp,ydnp,zdnp,lonnp,latnp)
        #debug
        xb, yb, zb = librotate.uv2xyzd(unp,vnp,lonnp,latnp)
        logging.debug(f"xyzd2uv x {np.max(np.abs(xdnp-xb))}")
        logging.debug(f"xyzd2uv y {np.max(np.abs(ydnp-yb))}")
        logging.debug(f"xyzd2uv z {np.max(np.abs(zdnp-zb))}")
        #missing value
        if(np.isnan(unp).sum() != 0 or np.isnan(vnp).sum() != 0):
            logging.warning("missing value exists in xyzd2uv")
            continue     
        udata = xr.DataArray(unp.reshape(1,nlat,nlon),\
                    [('time',pd.date_range(date,periods=1)),\
                    ('lat',latin),('lon',lonin)],\
                    attrs=attrs_u,name=var_sfc[0])
        vdata = xr.DataArray(vnp.reshape(1,nlat,nlon),\
                    [('time',pd.date_range(date,periods=1)),\
                    ('lat',latin),('lon',lonin)],\
                    attrs=attrs_v,name=var_sfc[1])
        logging.debug(udata)
        logging.debug(vdata)
        daout.append(udata)
        daout.append(vdata)
        
        #pl
        datain = xr.open_dataset(datadir / inplnc).sel(time=date)
        logging.debug(datain)
        lev = datain["level"].values
        nlev = len(lev)
        logging.debug(f"level: {lev}")
        uname = var_pl[0]
        vname = var_pl[1]
        #print(data[uname])
        u = datain[uname].values
        v = datain[vname].values 
        attrs_u = datain[uname].attrs
        attrs_v = datain[vname].attrs
        #missing values need to be searched
        if(np.isnan(u).sum() != 0 or np.isnan(v).sum() != 0):
            logging.warning("missing value exists in input")
            continue
        #   lon2d = lon.reshape(1,-1)
        #   lat2d = lat.reshape(-1,1)
        lon3d = lon[None, None, :]
        lat3d = lat[None, :, None]
        xd,yd,zd = librotate.uv2xyzd(u,v,lon3d,lat3d)
        #debug
        ud, vd = librotate.xyzd2uv(xd,yd,zd,lon3d,lat3d)
        logging.debug(f"uv2xyzd u {np.max(np.abs(u-ud))}")
        logging.debug(f"uv2xyzd v {np.max(np.abs(v-vd))}")
        xb, yb, zb = librotate.uv2xyzd(ud,vd,lon3d,lat3d)
        logging.debug(f"xyzd2uv x {np.max(np.abs(xd-xb))}")
        logging.debug(f"xyzd2uv y {np.max(np.abs(yd-yb))}")
        logging.debug(f"xyzd2uv z {np.max(np.abs(zd-zb))}")

        da_xd = xr.DataArray(xd.reshape(1,nlev,len(lat),len(lon)),\
                            [('time',pd.date_range(date,periods=1)),\
                            ('level',lev),\
                            ('latitude',lat),('longitude',lon)],\
                            name='xd')
        da_yd = xr.DataArray(yd.reshape(1,nlev,len(lat),len(lon)),\
                            [('time',pd.date_range(date,periods=1)),\
                            ('level',lev),    
                            ('latitude',lat),('longitude',lon)],\
                            name='yd')
        da_zd = xr.DataArray(zd.reshape(1,nlev,len(lat),len(lon)),\
                            [('time',pd.date_range(date,periods=1)),\
                            ('level',lev),
                            ('latitude',lat),('longitude',lon)],\
                            name='zd')
        
        xd_interp = da_xd.interp(longitude=newlon, latitude=newlat)
        yd_interp = da_yd.interp(longitude=newlon, latitude=newlat)
        zd_interp = da_zd.interp(longitude=newlon, latitude=newlat)
        #missing value
        #        dfx = xd_interp.to_pandas()
        #        dfy = yd_interp.to_pandas()
        #        dfz = zd_interp.to_pandas()
        dfx = xd_interp.values
        dfy = yd_interp.values
        dfz = zd_interp.values
        if(np.isnan(dfx).sum() != 0 or np.isnan(dfy).sum() != 0\
            or np.isnan(dfz).sum() != 0):
            logging.warning("missing value exists in interpolation")
            msk = np.isnan(dfx[0,:])
            mislon = newlon.values[msk]
            mislat = newlat.values[msk]
            for i in range(len(mislon)):
                logging.warning("missing at lon{} lat{}".format(mislon[i],mislat[i]))
            continue
        xdtc = xd_interp.values.reshape(1,nlev,nlat,nlon)
        ydtc = yd_interp.values.reshape(1,nlev,nlat,nlon)
        zdtc = zd_interp.values.reshape(1,nlev,nlat,nlon)
        xdnp,ydnp,zdnp = librotate.tc2np(lonc,latc,xdtc,ydtc,zdtc)
        #debug
        xb,yb,zb = librotate.np2tc(lonc,latc,xdnp,ydnp,zdnp)
        logging.debug(f"tc2np x {np.max(np.abs(xdtc-xb))}")
        logging.debug(f"tc2np y {np.max(np.abs(ydtc-yb))}")
        logging.debug(f"tc2np z {np.max(np.abs(zdtc-zb))}")
        #missing value
        if(np.isnan(xdnp).sum() != 0 or np.isnan(ydnp).sum() != 0 \
            or np.isnan(zdnp).sum() != 0):
            logging.warning("missing value exists in tc2np")
            continue
        #lonnp = lonin.reshape(1,-1)
        lonnp = lonin[None, None, None, :]
        latnp = latin[None, None, :, None]
            #unp,vnp = librotate.xyzd2uv(xdnp,ydnp,zdnp,lonnp)
        unp,vnp = librotate.xyzd2uv(xdnp,ydnp,zdnp,lonnp,latnp)
        #debug
        xb, yb, zb = librotate.uv2xyzd(unp,vnp,lonnp,latnp)
        logging.debug(f"xyzd2uv x {np.max(np.abs(xdnp-xb))}")
        logging.debug(f"xyzd2uv y {np.max(np.abs(ydnp-yb))}")
        logging.debug(f"xyzd2uv z {np.max(np.abs(zdnp-zb))}")
        #missing value
        if(np.isnan(unp).sum() != 0 or np.isnan(vnp).sum() != 0):
            logging.warning("missing value exists in xyzd2uv")
            continue        
        udata = xr.DataArray(unp,\
                    [('time',pd.date_range(date,periods=1)),\
                    ('level',lev),\
                    ('lat',latin),('lon',lonin)],\
                    attrs=attrs_u,name=uname)
        vdata = xr.DataArray(vnp,\
                    [('time',pd.date_range(date,periods=1)),\
                    ('level',lev),\
                    ('lat',latin),('lon',lonin)],\
                    attrs=attrs_v,name=vname)
        logging.debug(udata)
        logging.debug(vdata)
        daout.append(udata)
        daout.append(vdata)

        da_np.append(xr.merge(daout))
data_np = xr.merge(da_np)
logging.info(data_np)

data_np.to_netcdf(outdir / outnc, 'w')
