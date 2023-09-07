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
#vtime = datetime(2019, 10, 13, 3) # Hagibis changed to extra-tropical cyclone
dt = timedelta(hours = 12) # 1 days, 0:00:00
dc = 5
pc = 101000
latc = 40
sigma = 0.0
orig = "rjtd"
Cendict={"ecmf":"ecmwf","rjtd":"jma","kwbc":"ncep","egrr":"ukmo"}
center = Cendict[orig]
Memdict={"ecmwf":50,"jma":26,"ncep":20,"ukmo":17} #old
#MemList={"ecmwf":50,"jma":50,"ncep":30,"ukmo":17} #latest
Mem = Memdict[center]
#outdir = Path("../ecmwf")
#outdir.mkdir(0o755, True, True)
#for m in range(Mem):
t0 = datetime(2019, 10, 7, 0) # 2019-10-06 12:00:00
t0max = datetime(2019, 10, 12, 12) # 2019-10-07 12:00:00
while t0 <= t0max:
    init = t0.strftime("%m%d%H") # -> 2019100600
    year = t0.strftime("%Y")
    month = t0.strftime("%m")
    print(init)
    #outfile = orig+"/track" + year + init + "_mean.txt"
    outfile = center+"/track" + year + init + "_cntl.txt"
        #outfile = center+"/track" + year + init + "_" + str(m+1).zfill(2) + ".txt"
    #outfile = orig+"/track_anl.txt"
    track = open(outfile, "w")
    #infile = "/Volumes/dandelion/netcdf/tigge/"+year+"/"+center+"/"+init+"_mean.nc"
    infile = "/Volumes/dandelion/netcdf/tigge/"+year+"/"+month+"/"+center+"/glb_"+year+init+"_cntl.nc"
        #infile = "/Volumes/dandelion/netcdf/tigge/"+year+"/"+center+"/"+init+"_" + str(m+1).zfill(2) + ".nc"
    #infile = "../../netcdf/tigge/"+year+"/"+orig+"/anl.nc"
    nc = netCDF4.Dataset(infile, 'r')
    lon = nc.variables['longitude'][:]
    lat = nc.variables['latitude'][:]
    time = nc.variables['time']
    lonpre = None 
    latpre = None
    for i in range(len(time)):
        t = netCDF4.num2date(time[i], time.units) # 2019-10-06 12:00:00
        print(t)
        #if t > vtime: break
        if t > btime[-1]: break
        index_t0 = btime.index(t)
        lonb = bst[index_t0, 4]
        latb = bst[index_t0, 5]
        slpb = bst[index_t0, 6] #[Pa]
        #if i == 0:
        lon0 = lonb
        lat0 = latb
        slp0 = slpb
        print("best track center")
        print(f"{lonb:.3f}, {latb:.3f}, {slpb:.3f}")
        #slp = nc.variables['PRES_meansealevel'][i,]
        slp = nc.variables['msl'][i,]
        #if nc.variables['PRES_meansealevel'].units == 'hPa':
        if nc.variables['msl'].units == 'hPa':
            slp *= 1e2
        lonmin, latmin, slpmin = \
        grid.find_minimum(slp, lon, lat, lon0, lat0, slp0, lonpre, latpre, sigma)
        print("estimated center")
        print(f"{lonmin:.3f}, {latmin:.3f}, {slpmin:.3f}")
#        y0 = np.deg2rad(lat0)
#        y1 = np.deg2rad(latmin)
#        dx = np.deg2rad(lonmin - lon0)
#        d = np.arccos(np.sin(y0) * np.sin(y1) + np.cos(y0) * np.cos(y1) * np.cos(dx))
#        if np.rad2deg(d) > dc: break
        if slpmin > pc or latmin > latc : break
        lonpre = lonmin
        latpre = latmin
        print("{} {} {} {} {} {} {}".format(t.year, t.month, t.day, t.hour, lonmin, latmin, slpmin), file = track)
    nc.close()
    t0 += dt



