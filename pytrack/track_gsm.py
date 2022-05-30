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

dt = timedelta(hours = 12) # 1 days, 0:00:00
vkey = "Z300"
dc = 5
if vkey == "PSEA":
    pc = 1010
    plim = 20.0
else:
    pc = 10000
    plim = 1000.0
    plev = float(vkey[1:])
latc = 40
sigma = 0.0
#outdir = Path("../ecmwf")
#outdir.mkdir(0o755, True, True)
suffix = "_mgdsst"
t0 = datetime(2019, 10, 9, 0) # 2019-10-06 12:00:00
t0max = datetime(2019, 10, 9, 12) # 2019-10-07 12:00:00
while t0 <= t0max:
    init = t0.strftime("%Y%m%d%H") # -> 2019100600
    year = t0.strftime("%Y")
    month = t0.strftime("%m")
    day = t0.strftime("%d")
    hour = t0.strftime("%H")
    print(init)
    if vkey == "PSEA":
        outfile = "./track" + init + "_gsm_tl959"+suffix+".txt"
        infile = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+suffix+"/fcst_surf_" + init + ".nc"
        #infile = "/Volumes/dandelion/GSMJob/Jobwk_Tl319L100_sstclim+mysub/fcst_p_" + init + ".nc"
    else:
        outfile = "./track" + vkey[1:] + "_" + init + "_gsm_tl959"+suffix+".txt"
        infile = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+suffix+"/fcst_p_asia_" + init + ".nc"
    track = open(outfile, "w")
    nc = netCDF4.Dataset(infile, 'r')
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    time = nc.variables['time']
    for i in range(len(time)):
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
        print(lon0, lat0)
        if vkey == "PSEA":
            slp = nc.variables[vkey][i,]
        else:
            level = nc.variables['level'][:]
            ind_p = np.argmin(np.abs(level - plev))
            slp = nc.variables[vkey[0]][i,ind_p,:,:]
        lonmin, latmin, slpmin = grid.find_minimum(slp, lon, lat, lon0, lat0, sigma)
        print(slpmin)
        if slpmin > pc and latmin > latc : break
        # if track is far away from the previous center, ignore
        d = 0.0
        p = 0.0
        if i > 0:
            y0 = np.deg2rad(latpre)
            y1 = np.deg2rad(latmin)
            dx = np.deg2rad(lonmin - lonpre)
            d = np.arccos(np.sin(y0) * np.sin(y1) + np.cos(y0) * np.cos(y1) * np.cos(dx))
            p = np.abs(slpmin - slppre)
        print(d,p)
        if np.rad2deg(d) > dc: continue
        if p > plim: continue
        lonpre = lonmin
        latpre = latmin
        slppre = slpmin
        #if t.hour % 6 == 0:
        print("{} {} {} {} {} {} {}".format(t.year, t.month, t.day, t.hour, lonmin, latmin, slpmin), file = track)
    nc.close()
    t0 += dt



