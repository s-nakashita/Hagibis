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

dt = timedelta(hours = 24) # 1 days, 0:00:00
dc = 5
pc = 101000
latc = 40
sigma = 0.0
orig = "jma"
Mem = 26
#outdir = Path("../ecmwf")
#outdir.mkdir(0o755, True, True)
#for m in range(Mem):
t0 = datetime(2019, 10, 9, 0) # 2019-10-06 12:00:00
t0max = datetime(2019, 10, 13, 0) # 2019-10-07 12:00:00
outfile = "track_era5.txt"
        #outfile = orig+"/track" + year + init + "_" + str(m+1).zfill(2) + ".txt"
    #outfile = orig+"/track_anl.txt"
track = open(outfile, "w")
while t0 <= t0max:
    init = t0.strftime("%m%d%H") # -> 2019100600
    year = t0.strftime("%Y")
    mm   = t0.strftime("%m")
    dd   = t0.strftime("%d")
    print(init)
    infile = "/Volumes/dandelion/netcdf/era5/"+year+"/"+mm+"/"+dd+"/slev_"+dd+".nc"
        #infile = "../../netcdf/tigge/"+year+"/"+orig+"/"+init+"_" + str(m+1).zfill(2) + ".nc"
    nc = netCDF4.Dataset(infile, 'r')
    lon = nc.variables['longitude'][:]
    lat = nc.variables['latitude'][:]
    time = nc.variables['time']
    #for i in range(len(time)):
    for i in range(8):
        t = netCDF4.num2date(time[i], time.units) # 2019-10-06 12:00:00
        print(t)
        if t > btime[-1]: break
        try:
            index_t0 = btime.index(t)
        except ValueError:
            # use average between prior and post values
            t_pre = t - timedelta(hours = 3)
            t_pos = t + timedelta(hours = 3)
            index_tpre = btime.index(t_pre)
            index_tpos = btime.index(t_pos)
            lon1 = bst[index_tpre, 4]
            lat1 = bst[index_tpre, 5]
            lon2 = bst[index_tpos, 4]
            lat2 = bst[index_tpos, 5]
            lon0 = (lon1 + lon2) / 2
            lat0 = (lat1 + lat2) / 2
            #continue
        else:
            lon0 = bst[index_t0, 4]
            lat0 = bst[index_t0, 5]
        print("best track center")
        print(lon0, lat0)
        slp = nc.variables['msl'][i,]
        lonmin, latmin, slpmin = grid.find_minimum(slp, lon, lat, lon0, lat0, sigma)
        print("estimated center")
        print(lonmin, latmin)
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



