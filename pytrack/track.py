import netCDF4
import grid
import sys
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np

# year mon day hour lon   lat  slp
# 2019 10  4   18   164.4 15.7 100800
bst = np.loadtxt("bst_hagibis.txt")
btime = []
for t in bst[:, (0, 1, 2, 3)].astype(np.int32): # [2019 10 4 18]
    btime.append(datetime(*t)) # 2019-10-04 18:00:00

t0 = datetime(2019, 10, 9, 0) # 2019-10-06 12:00:00
t0max = datetime(2019, 10, 9, 12) # 2019-10-07 12:00:00
dt = timedelta(hours = 12) # 1 days, 0:00:00
dc = 5
pc = 101000
latc = 40
sigma = 0.0
orig = "jma"
#outdir = Path("../ecmwf")
#outdir.mkdir(0o755, True, True)

while t0 <= t0max:
    init = t0.strftime("%m%d%H") # -> 2019100600
    year = t0.strftime("%Y")
    print(init)
    outfile = orig+"/track" + year + init + "_c.txt"
    track = open(outfile, "w")
    infile = "../../netcdf/tigge/"+year+"/"+orig+"/"+init+"_ctrl.nc"
    nc = netCDF4.Dataset(infile, 'r')
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    time = nc.variables['time']
    for i in range(len(time)):
        t = netCDF4.num2date(time[i], time.units) # 2019-10-06 12:00:00
        print(t)
        if t > btime[-1]: break
        index_t0 = btime.index(t)
        lon0 = bst[index_t0, 4]
        lat0 = bst[index_t0, 5]
        print(lon0, lat0)
        slp = nc.variables['PRES_meansealevel'][i,]
        lonmin, latmin, slpmin = grid.find_minimum(slp, lon, lat, lon0, lat0, sigma)
        print(slpmin)
#        y0 = np.deg2rad(lat0)
#        y1 = np.deg2rad(latmin)
#        dx = np.deg2rad(lonmin - lon0)
#        d = np.arccos(np.sin(y0) * np.sin(y1) + np.cos(y0) * np.cos(y1) * np.cos(dx))
#        if np.rad2deg(d) > dc: break
        if slpmin > pc and latmin > latc : break
        lon0 = lonmin
        lat0 = latmin
        print("{} {} {} {} {} {} {}".format(t.year, t.month, t.day, t.hour, lonmin, latmin, slpmin), file = track)
    nc.close()
    t0 += dt



