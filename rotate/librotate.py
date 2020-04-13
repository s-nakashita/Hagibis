import sys
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import xarray as xr

pi = np.pi
dtheta = pi/180.0
#lon,lat is degree

def lonlat2xyz(lon,lat):
    x = np.cos(lon*dtheta) * np.cos(lat*dtheta)
    y = np.sin(lon*dtheta) * np.cos(lat*dtheta)
    z = np.sin(lat*dtheta)
    return x,y,z

def xyz2lonlat(x,y,z):
    if(x != 0.0):
        lon = np.arctan2(y,x)
        if(lon < 0.0):
            lon += 2.0*pi
    else:
        lon = 0.0
    lon /= dtheta
    lat = np.arcsin(z)/dtheta
    return lon,lat

def uv2xyzd(u,v,lon,lat):
    xd = -u*np.sin(lon*dtheta) - v*np.cos(lon*dtheta)*np.sin(lat*dtheta)
    yd = u*np.cos(lon*dtheta) - v*np.sin(lon*dtheta)*np.sin(lat*dtheta)
    zd = v*np.cos(lat*dtheta)
    return xd,yd,zd

def xyzd2uv(xd,yd,zd,lon):
    u = -xd*np.sin(lon*dtheta) + yd*np.cos(lon*dtheta)
    v = sgn(zd) * np.sqrt((xd*np.cos(lon*dtheta) + yd*np.sin(lon*dtheta))**2 + zd**2)
    return u,v

def generate_points(nlon,nlat,dlat):
    dlon = 360.0/nlon
    lon = dlon * np.arange(nlon)
    lat = 90.0 - dlat * np.arange(nlat)
    return lon,lat

#Ax
def np2tc(lonc,latc,x,y,z):
    xx = np.cos(lonc*dtheta)*np.sin(latc*dtheta)*x\
        - np.sin(lonc*dtheta)*y\
        + np.cos(lonc*dtheta)*np.cos(latc*dtheta)*z
    yy = np.sin(lonc*dtheta)*np.sin(latc*dtheta)*x\
        + np.cos(lonc*dtheta)*y\
        + np.sin(lonc*dtheta)*np.cos(latc*dtheta)*z
    zz = -np.cos(latc*dtheta)*x + np.sin(latc*dtheta)*z
    return xx,yy,zz

#ATx
def tc2np(lonc,latc,x,y,z):
    xx = np.cos(lonc*dtheta)*np.sin(latc*dtheta)*x\
        + np.sin(lonc*dtheta)*np.sin(latc*dtheta)*y\
        - np.cos(latc*dtheta)*z
    yy = -np.sin(lonc*dtheta)*x + np.cos(lonc*dtheta)*y
    zz = np.cos(lonc*dtheta)*np.cos(latc*dtheta)*x\
        + np.sin(lonc*dtheta)*np.cos(latc*dtheta)*y\
        + np.sin(latc*dtheta)*z
    return xx,yy,zz

def rotate_lonlat(lonc,latc,lon,lat):
    nlon = lon.size
    nlat = lat.size
    lonout = np.zeros(nlon*nlat)
    latout = np.zeros(nlon*nlat)
    ij = 0
    for j in range(nlat):
        for i in range(nlon):
            x,y,z = lonlat2xyz(lon[i],lat[j])
            xx,yy,zz = np2tc(lonc,latc,x,y,z)
            lonout[ij],latout[ij] = xyz2lonlat(xx,yy,zz)
            ij += 1
    return lonout,latout

def sgn(x):
    return (x>0).astype(np.int) - (x<0).astype(np.int)
