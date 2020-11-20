import sys
import numpy as np
from pathlib import Path
from datetime import datetime,timedelta
import xarray as xr
import pandas as pd
import librotate

#Usage echo yyyymmddhh datadir trackf nlon nlat latmax | python rotate_uv.py
param = sys.stdin.readline().strip("\n").split(" ")
yyyymmddhh = param[0]
ddirname   = param[1]
trackname  = param[2]
nlon       = int(param[3])
nlat       = int(param[4])
latmax     = float(param[5]) #degree
dlat       = latmax/(nlat-1) #degree

datadir = Path(ddirname)
outdir  = Path('./')
mmddhh  = yyyymmddhh[4:]
innc    = mmddhh + '_mean.nc'
outnc   = 'np_ve_' + yyyymmddhh + '_mean.nc'
#innc    = 'anl.nc'
#outnc   = 'np_ve_anl.nc'
trackf  = Path(trackname)

var_sfc = ['UGRD_10maboveground','VGRD_10maboveground']
var_pl  = ['UGRD','VGRD']

lonin,latin = librotate.generate_points(nlon,nlat,dlat)

da_np = []
with trackf.open() as track:
    for l in track:
        l_split = l.split()
        print(l_split)
        year     = l_split[0]
        month    = l_split[1]
        day      = l_split[2]
        hour     = l_split[3]
        date_str = year + '/' + month + '/' + day + ' ' + hour + ':00'
        date     = datetime.strptime(date_str, '%Y/%m/%d %H:%M')
        print(date)
        lonc     = float(l_split[4])
        latc     = float(l_split[5])
        print(lonc,latc)
        
        lonout,latout = librotate.rotate_lonlat(lonc,latc,lonin,latin)
        print(lonout.max(),lonout.min(),latout.max(),latout.min())
        
        data = xr.open_dataset(datadir / innc).sel(time=date)
        print(data)
        
        lon = data["lon"].values
        lat = data["lat"].values
        lev = data["level"].values
        nlev = len(lev)
        newlon = xr.DataArray(lonout,dims='np_lonlat')
        newlat = xr.DataArray(latout,dims='np_lonlat')
        daout = []
        
        #sfc
        u = data[var_sfc[0]].values
        v = data[var_sfc[1]].values 
        attrs_u = data[var_sfc[0]].attrs
        attrs_v = data[var_sfc[1]].attrs
        #missing values need to be searched
        dfu = data[var_sfc[0]].to_pandas()
        dfv = data[var_sfc[1]].to_pandas()
        if(dfu.isnull().values.sum() != 0 or dfv.isnull().values.sum() != 0):
            print("missing value exists in input")
        else:
            lon2d = lon.reshape(1,-1)
            lat2d = lat.reshape(-1,1)
            xd,yd,zd = librotate.uv2xyzd(u,v,lon2d,lat2d)
        
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
            dfx = xd_interp.to_pandas()
            dfy = yd_interp.to_pandas()
            dfz = zd_interp.to_pandas()
            if(dfx.isnull().values.sum() != 0 or dfy.isnull().values.sum() != 0\
               or dfz.isnull().values.sum() != 0):
                print("missing value exists in interpolation")
            else:
                xdtc = xd_interp.values.reshape(nlat,nlon)
                ydtc = yd_interp.values.reshape(nlat,nlon)
                zdtc = zd_interp.values.reshape(nlat,nlon)
                xdnp,ydnp,zdnp = librotate.tc2np(lonc,latc,xdtc,ydtc,zdtc)
                #missing value
                if(np.isnan(xdnp).sum() != 0 or np.isnan(ydnp).sum() != 0 or \
                   np.isnan(zdnp).sum() != 0):
                    print("missing value exists in tc2np")
                else:
                    lonnp = lonin.reshape(1,-1)
                    unp,vnp = librotate.xyzd2uv(xdnp,ydnp,zdnp,lonnp)
                    #missing value
                    if(np.isnan(unp).sum() != 0 or np.isnan(vnp).sum() != 0):
                        print("missing value exists in xyzd2uv")
                    else:        
                        udata = xr.DataArray(unp.reshape(1,nlat,nlon),\
                             [('time',pd.date_range(date,periods=1)),\
                              ('latitude',latin),('longitude',lonin)],\
                             attrs=attrs_u,name=var_sfc[0])
                        vdata = xr.DataArray(vnp.reshape(1,nlat,nlon),\
                             [('time',pd.date_range(date,periods=1)),\
                              ('latitude',latin),('longitude',lonin)],\
                             attrs=attrs_v,name=var_sfc[1])
                        print(udata)
                        print(vdata)
                        daout.append(udata)
                        daout.append(vdata)
        
        #pl
        #for lev in ['200','250','300','500','850']:
        #uname = var_pl[0].replace('000',lev)
        #vname = var_pl[1].replace('000',lev)
        uname = var_pl[0]
        vname = var_pl[1]
        #print(data[uname])
        level = []
        for l in range(nlev):
            u = data[uname].values[l]
            v = data[vname].values[l] 
            attrs_u = data[uname].attrs
            attrs_v = data[vname].attrs
            #missing values need to be searched
            dfu = data[uname][l,:].to_pandas()
            dfv = data[vname][l,:].to_pandas()
            if(dfu.isnull().values.sum() != 0 or dfv.isnull().values.sum() != 0):
                print("missing value exists in input")
            else:
                lon2d = lon.reshape(1,-1)
                lat2d = lat.reshape(-1,1)
                xd,yd,zd = librotate.uv2xyzd(u,v,lon2d,lat2d)
        
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
                dfx = xd_interp.to_pandas()
                dfy = yd_interp.to_pandas()
                dfz = zd_interp.to_pandas()
                if(dfx.isnull().values.sum() != 0 or dfy.isnull().values.sum() != 0\
                   or dfz.isnull().values.sum() != 0):
                    print("missing value exists in interpolation")
                else:
                    xdtc = xd_interp.values.reshape(nlat,nlon)
                    ydtc = yd_interp.values.reshape(nlat,nlon)
                    zdtc = zd_interp.values.reshape(nlat,nlon)
                    xdnp,ydnp,zdnp = librotate.tc2np(lonc,latc,xdtc,ydtc,zdtc)
                    #missing value
                    if(np.isnan(xdnp).sum() != 0 or np.isnan(ydnp).sum() != 0 \
                       or np.isnan(zdnp).sum() != 0):
                        print("missing value exists in tc2np")
                    else:
                        lonnp = lonin.reshape(1,-1)
                        unp,vnp = librotate.xyzd2uv(xdnp,ydnp,zdnp,lonnp)
                        #missing value
                        if(np.isnan(unp).sum() != 0 or np.isnan(vnp).sum() != 0):
                            print("missing value exists in xyzd2uv")
                        else:        
                            udata = xr.DataArray(unp.reshape(1,1,nlat,nlon),\
                                    [('time',pd.date_range(date,periods=1)),\
                                     ('level',np.array([lev[l]])),\
                                     ('latitude',latin),('longitude',lonin)],\
                                                 attrs=attrs_u,name=uname)
                            vdata = xr.DataArray(vnp.reshape(1,1,nlat,nlon),\
                                    [('time',pd.date_range(date,periods=1)),\
                                     ('level',np.array([lev[l]])),\
                                     ('latitude',latin),('longitude',lonin)],\
                                                 attrs=attrs_v,name=vname)
                        print(udata)
                        print(vdata)
                        daout.append(udata)
                        daout.append(vdata)
                        
        da_np.append(xr.merge(daout))
data_np = xr.merge(da_np)
print(data_np)

data_np.to_netcdf(outdir / outnc, 'w')
