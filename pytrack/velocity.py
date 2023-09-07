import sys
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
from numpy import pi, sin, cos, tan, arcsin, arccos, arctan

deg2rad = pi / 180.0
rad2deg = 180.0 / pi
re = 6.371e6

def speed(lon1, lon2, lat1, lat2, dt):
    l1 = lon1 * deg2rad
    l2 = lon2 * deg2rad
    p1 = lat1 * deg2rad
    p2 = lat2 * deg2rad

    return re * arccos(sin(p1)*sin(p2) + cos(p1)*cos(p2)*cos(l1 - l2)) / dt

def direction(lon1, lon2, lat1, lat2):
    l1 = lon1 * deg2rad
    l2 = lon2 * deg2rad
    p1 = lat1 * deg2rad
    p2 = lat2 * deg2rad

    vl = arccos(sin(p1)**2 + cos(p1)**2*cos(l1 - l2))
    vp = arccos(sin(p1)*sin(p2) + cos(p1)*cos(p2))
    if l1 < l2 and p1 < p2:
        theta = arctan(vl / vp)
    elif l1 > l2 and p1 < p2:
        theta = -arctan(vl / vp)
    elif l1 > l2 and p1 > p2:
        theta = -0.5*pi + arctan(vl / vp)
    else:
        theta = 0.5*pi - arctan(vl / vp)
    return theta * rad2deg
# year mon day hour lon   lat  slp
# 2019 10  4   18   164.4 15.7 100800
outbst = False
if outbst:
    bst = np.loadtxt("bst_hagibis.txt")
    btime = []
    ind = []
    i = -1
    for t in bst[:, (0, 1, 2, 3)].astype(np.int32): # [2019 10 4 18]
        i += 1
        if not t[3] % 6 == 0:
            print(t[3])
            continue 
        btime.append(datetime(*t)) # 2019-10-04 18:00:00
        ind.append(i)
    print(bst[(ind),3])
    outfile = "bst_velocity.txt"
    velocity = open(outfile, "w")
    for l in range(len(btime)):
        ind1 = max(0,l-1)
        ind2 = min(len(btime)-1,l+1)
        lon1 = bst[ind[ind1], 4]
        lon2 = bst[ind[ind2], 4]
        lat1 = bst[ind[ind1], 5]
        lat2 = bst[ind[ind2], 5]
        dt = (btime[ind2] - btime[ind1]).total_seconds()
        v = speed(lon1, lon2, lat1, lat2, dt)
        d = direction(lon1, lon2, lat1, lat2)
        print("{} {} {} {} {} {}".format(int(btime[l].year), int(btime[l].month), \
        int(btime[l].day), int(btime[l].hour), v, d),\
         file=velocity)

origin = ["ecmf", "rjtd", "kwbc", "egrr"]
center = {"ecmf":"ecmwf","rjtd":"jma","kwbc":"ncep","egrr":"ukmo"}
dt0 = timedelta(hours = 12)
for orig in ["rjtd"]:
    t0 = datetime(2019, 10, 10, 12)
    t0max = datetime(2019, 10, 10, 12)
    while t0 <= t0max:
        init = t0.strftime("%Y%m%d%H")
        yyyy = t0.strftime("%Y")
        tfile = center[orig] + "/track" + init + "_cntl.txt"
        ofile = center[orig] + "/velocity" + init + "_cntl.txt"
        track = np.loadtxt(tfile)
        velocity = open(ofile, "w")
        time = []
        ind = []
        i = -1
        for t in track[:, (0, 1, 2, 3)].astype(np.int32): # [2019 10 4 18]
            i += 1
            if not t[3] % 6 == 0:
                print(t[3])
                continue 
            time.append(datetime(*t)) # 2019-10-04 18:00:00
            ind.append(i)
        print(track[(ind),3])
        for l in range(len(time)):
            ind1 = max(0,l-1)
            ind2 = min(len(time)-1,l+1)
            lon1 = track[ind[ind1], 4]
            lon2 = track[ind[ind2], 4]
            lat1 = track[ind[ind1], 5]
            lat2 = track[ind[ind2], 5]
            dt = (time[ind2] - time[ind1]).total_seconds()
            v = speed(lon1, lon2, lat1, lat2, dt)
            d = direction(lon1, lon2, lat1, lat2)
            vwe = v*sin(d*deg2rad)
            vns = v*cos(d*deg2rad)
            print("{} {} {} {} {} {} {} {} {} {}".format(\
                int(time[l].year), int(time[l].month), \
                int(time[l].day), int(time[l].hour), \
                track[ind[l], 4], track[ind[l], 5],\
                v, d, vwe, vns),\
                file=velocity)
        t0 += dt0