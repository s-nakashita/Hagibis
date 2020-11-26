import sys
import numpy as np
from pathlib import Path
from datetime import datetime,timedelta
import logging
import xarray as xr
import pandas as pd
import librotate

logging.basicConfig(level=logging.INFO)

#Usage python rotate_uv.py yyyymmddhh datadir trackf dcolat
#param = sys.stdin.readline().strip("\n").split(" ")
yyyymmddhh = sys.argv[1]
ddirname   = sys.argv[2]
trackname  = sys.argv[3]
dcolat     = float(sys.argv[4]) #degree
logging.info(dcolat)

datadir = Path(ddirname)
outdir  = Path('./')
mmddhh  = yyyymmddhh[4:]
innc    = 'np_glb_' + yyyymmddhh + '_mean.nc'
outnc   = 'inv1d_' + yyyymmddhh + '_mean.nc'
#innc    = 'anl.nc'
#outnc   = 'np_ve_anl.nc'
trackf  = Path(trackname)
logging.debug(trackf)

var_scl = ['msl']
var_pl  = ['u','v']

dain = xr.open_dataset(datadir / innc) # Dataset
lon_attrs = dain["lon"].attrs
lat_attrs = dain["lat"].attrs
# add cyclic
dabc = dain.sel(lon=0.0)
dabc["lon"] = 360.0
logging.debug(dabc)
da = xr.concat([dain, dabc], dim="lon")
logging.debug(da)
lonin = da["lon"].values
latin = da["lat"].values
level = da["level"].values
nlev = level.size
logging.debug(f"level: {level}")
time = da["time"]
tmax = time[-1]
logging.info(time)

# mask area
lat_msk = da["lat"].where(da.lat >= 90.0 - dcolat, drop=True).values.tolist()
logging.info(lat_msk)
#quit()
# global grid (every 1 degree)
lon = np.arange(360)
lat = -90.0 + np.arange(181)
nlon = lon.size
nlat = lat.size
logging.debug(lon)
logging.debug(lat)

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
        logging.debug(time.isin(pd.date_range(date,periods=1)).values.sum())
        #exit()
        if time.isin(pd.date_range(date,periods=1)).values.sum() == 0:
            break
        lonc     = float(l_split[4])
        latc     = float(l_split[5])
        logging.info(f"TC center ({lonc},{latc})")

        data = da.sel(time=date)
        logging.debug(data["lat"])

        daout = []
        
        # create interpolated grid
        msk = librotate.mask_lonlat(60.0, lonc, latc, lon, lat)
        nmsk = np.sum(msk.astype(np.int))
        logging.info(f"output grid size: {nmsk}")
        ndim = xr.DataArray(int(nmsk),[('time',pd.date_range(date,periods=1))],\
            name="GridSize")
        daout.append(ndim)
        lontc = np.zeros(nmsk)
        lattc = np.zeros(nmsk)
        k = 0
        for j in range(nlat):
            for i in range(nlon):
                if msk[j, i]:
                    lontc[k] = lon[i]
                    lattc[k] = lat[j]
                    k += 1
        logging.debug(\
            f"output coordinate \n \
            lon:min{lontc.min():.3f}, max{lontc.max():.3f} \n \
            lat:min{lattc.min():.3f}, max{lattc.max():.3f}")
        lonnp,latnp = librotate.rotate_lonlat1d(lonc,latc,lontc,lattc,-1)
        logging.debug(\
            f"interpolate coordinate \n \
            lon:min{lonnp.min():.3f}, max{lonnp.max():.3f} \n \
            lat:min{latnp.min():.3f}, max{latnp.max():.3f}")
        
        newlon = xr.DataArray(lonnp,dims='np_lonlat')
        newlat = xr.DataArray(latnp,dims='np_lonlat')
        
        #exit()
        
        # mask data
        # msl
        msl = data["msl"].values
        msl_mask = data["msl"].sel(lat=lat_msk)
        logging.debug(msl.shape)
        logging.debug(msl_mask.values.shape)
        msl_ = np.mean(msl_mask.values, axis=1) - np.mean(msl_mask.values[-1,:]) 
        logging.debug(msl_)
        msl[:len(lat_msk),:] = msl[:len(lat_msk),:] - msl_[:, None]
        logging.debug(msl)
        # u
        u = data["u"].values
        u_mask = data["u"].sel(lat=lat_msk)
        logging.debug(u.shape)
        logging.debug(u_mask.values.shape)
        u_ = np.mean(u_mask.values, axis=2) 
        logging.debug(u_)
        u[:, :len(lat_msk), :] = u[:, :len(lat_msk), :] - u_[:, :, None]
        logging.debug(u)
        # v
        v = data["v"].values
        v_mask = data["v"].sel(lat=lat_msk)
        logging.debug(v.shape)
        logging.debug(v_mask.values.shape)
        v_ = np.mean(v_mask.values, axis=2) 
        logging.debug(v_)
        v[:, :len(lat_msk), :] = v[:, :len(lat_msk), :] - v_[:, :, None]
        logging.debug(v)
        
        # rotation
        # msl
        attrs = data["msl"].attrs
        da_interp = data["msl"].interp(lon=newlon, lat=newlat)
        # missing value
        val = da_interp.values
        logging.debug(np.isnan(val).sum() )
        if np.isnan(val).sum() != 0:
            logging.warning("#{} missing value exists in interpolation"\
            .format(np.isnan(val).sum()))
            msk = np.isnan(val)
            mislon = newlon[msk]
            mislat = newlat[msk]
            for i in range(len(mislon)):
                logging.warning("missing at lon{} lat{}".format(mislon[i].values,mislat[i].values))
            continue
        lonout = xr.DataArray(lontc.reshape(1,-1),\
            [('time',pd.date_range(date,periods=1)),\
            ('lonlat1d',np.arange(len(lontc)))],\
                attrs=lon_attrs,name='lon1d')
        latout = xr.DataArray(lattc.reshape(1,-1),\
            [('time',pd.date_range(date,periods=1)),\
            ('lonlat1d',np.arange(len(lattc)))],\
                attrs=lat_attrs,name='lat1d')
        logging.debug(lonout)
        logging.debug(latout)
        daout.append(lonout)
        daout.append(latout)
        logging.debug(val.shape)
        #exit()
        #val = np.zeros((1,nlat,nlon))
        #k = 0
        #for j in range(nlat):
        #    for i in range(nlon):
        #        if msk[j, i]:
        #            val[0, j, i] = da_interp.values[k]
        #            k += 1
        #        else:
        #            val[0, j, i] = -999
        #logging.debug(val.size - np.isnan(val).sum())
        #exit()
        dataout = xr.DataArray(val.reshape(1,-1),\
            [('time',pd.date_range(date,periods=1)),\
            ('lonlat1d',np.arange(len(lattc)))],\
            attrs=attrs,name='msl')
        logging.debug(dataout)
        daout.append(dataout)
        #exit()
        
        #u, v
        attrs_u = data['u'].attrs
        attrs_v = data['v'].attrs
        #missing values need to be searched
        if(np.isnan(u).sum() != 0 or np.isnan(v).sum() != 0):
            logging.warning("missing value exists in input")
            continue
        lon3d = lonin[None, None, :]
        lat3d = latin[None, :, None]
        xd,yd,zd = librotate.uv2xyzd(u,v,lon3d,lat3d)
        #debug
        ud, vd = librotate.xyzd2uv(xd,yd,zd,lon3d,lat3d)
        logging.debug(f"uv2xyzd u {np.max(np.abs(u-ud))}")
        logging.debug(f"uv2xyzd v {np.max(np.abs(v-vd))}")
        xb, yb, zb = librotate.uv2xyzd(ud,vd,lon3d,lat3d)
        logging.debug(f"xyzd2uv x {np.max(np.abs(xd-xb))}")
        logging.debug(f"xyzd2uv y {np.max(np.abs(yd-yb))}")
        logging.debug(f"xyzd2uv z {np.max(np.abs(zd-zb))}")

        da_xd = xr.DataArray(xd[None,:],\
                            [('time',pd.date_range(date,periods=1)),\
                            ('level',level),\
                            ('latitude',latin),('longitude',lonin)],\
                            name='xd')
        da_yd = xr.DataArray(yd[None,:],\
                            [('time',pd.date_range(date,periods=1)),\
                            ('level',level),    
                            ('latitude',latin),('longitude',lonin)],\
                            name='yd')
        da_zd = xr.DataArray(zd[None,:],\
                            [('time',pd.date_range(date,periods=1)),\
                            ('level',level),
                            ('latitude',latin),('longitude',lonin)],\
                            name='zd')
        logging.debug(da_xd)
        #exit()
        xd_interp = da_xd.interp(longitude=newlon, latitude=newlat)
        yd_interp = da_yd.interp(longitude=newlon, latitude=newlat)
        zd_interp = da_zd.interp(longitude=newlon, latitude=newlat)
        #missing value
        #    dfx = xd_interp.to_pandas()
        #    dfy = yd_interp.to_pandas()
        #    dfz = zd_interp.to_pandas()
        dfx = np.isnan(xd_interp).astype(np.int)
        dfy = np.isnan(yd_interp).astype(np.int)
        dfz = np.isnan(zd_interp).astype(np.int)
        #    if(dfx.isnull().values.sum() != 0 or dfy.isnull().values.sum() != 0\
        #        or dfz.isnull().values.sum() != 0):
        if(np.sum(dfx) != 0 or np.sum(dfy) != 0 or np.sum(dfz) != 0):
            logging.warning("missing value exists in interpolation")
            continue
        logging.debug(xd_interp.values.shape)
        #exit()
        xdnp = xd_interp.values
        ydnp = yd_interp.values
        zdnp = zd_interp.values
        #xdnp = np.zeros((1,nlev,nlat,nlon))
        #ydnp = np.zeros_like(xdnp)
        #zdnp = np.zeros_like(xdnp)
        #k = 0
        #for j in range(nlat):
        #    for i in range(nlon):
        #        if msk[j, i]:
        #            xdnp[:, :, j, i] = xd_interp.values[:, :, k]
        #            ydnp[:, :, j, i] = yd_interp.values[:, :, k]
        #            zdnp[:, :, j, i] = zd_interp.values[:, :, k]
        #            k += 1
        #        else:
        #            xdnp[:, :, j, i] = 0.0
        #            ydnp[:, :, j, i] = 0.0
        #            zdnp[:, :, j, i] = 0.0
        #logging.debug(xdnp.size - np.isnan(xdnp).sum())
        #exit()
        xdtc,ydtc,zdtc = librotate.np2tc(lonc,latc,xdnp,ydnp,zdnp)
        #debug
        xb,yb,zb = librotate.tc2np(lonc,latc,xdtc,ydtc,zdtc)
        logging.debug(f"np2tc x {np.max(np.abs(xdnp-xb))}")
        logging.debug(f"np2tc y {np.max(np.abs(ydnp-yb))}")
        logging.debug(f"np2tc z {np.max(np.abs(zdnp-zb))}")
        #exit()
        #missing value
        if(np.isnan(xdtc).sum() != 0 or np.isnan(ydtc).sum() != 0 \
            or np.isnan(zdtc).sum() != 0):
            logging.warning("missing value exists in np2tc")
            continue
        #lonnp = lonin.reshape(1,-1)
        #lontc = lon[None, None, None, :]
        #lattc = lat[None, None, :, None]
            #unp,vnp = librotate.xyzd2uv(xdnp,ydnp,zdnp,lonnp)
        utc,vtc = librotate.xyzd2uv(xdtc,ydtc,zdtc,lontc[None, None, :],lattc[None, None, :])
        #debug
        xb, yb, zb = librotate.uv2xyzd(utc,vtc,lontc,lattc)
        logging.debug(f"xyzd2uv x {np.max(np.abs(xdtc-xb))}")
        logging.debug(f"xyzd2uv y {np.max(np.abs(ydtc-yb))}")
        logging.debug(f"xyzd2uv z {np.max(np.abs(zdtc-zb))}")
        #missing value
        if(np.isnan(utc).sum() != 0 or np.isnan(vtc).sum() != 0):
            logging.warning("missing value exists in xyzd2uv")
            continue   
        udata = xr.DataArray(utc, #.reshape(1,1,nlat,nlon),\
                    [('time',pd.date_range(date,periods=1)),\
                    ('level',level),\
                    ('lonlat1d',np.arange(len(lattc)))],\
                    attrs=attrs_u,name='u')
        vdata = xr.DataArray(vtc, #.reshape(1,1,nlat,nlon),\
                    [('time',pd.date_range(date,periods=1)),\
                    ('level',level),\
                    ('lonlat1d',np.arange(len(lattc)))],\
                    attrs=attrs_v,name='v')
        logging.debug(udata)
        logging.debug(vdata)
        daout.append(udata)
        daout.append(vdata)

        da_np.append(xr.merge(daout))
data_np = xr.merge(da_np)
logging.info(data_np)

data_np.to_netcdf(outdir / outnc, 'w')
