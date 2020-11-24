import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

bst = np.loadtxt("bst_velocity.txt")
btime = []
for t in bst[:, (0, 1, 2, 3)].astype(np.int32):
    btime.append(datetime(*t))
v_bst = bst[:, 4]
d_bst = bst[:, 5]
index_tmax = btime.index(datetime(2019, 10, 13, 6))
origin = ["ecmf", "rjtd", "kwbc", "egrr"]
center = {"ecmf":"ecmwf", "rjtd":"jma", "kwbc":"ncep", "egrr":"ukmo"}
mcolor = {"ecmf":"blue", "rjtd":"red", "kwbc":"green", "egrr":"lightgreen"}
t0 = datetime(2019, 10, 8, 0) # 2019-10-06 12:00:00
t0max = datetime(2019, 10, 10, 12) # 2019-10-07 12:00:00
dt0 = timedelta(hours = 12)
tmax = datetime(2019, 10, 12, 12)
while t0 <= t0max:
    init = t0.strftime("%Y%m%d%H") # -> 2019100600
    year = t0.strftime("%Y")
    print(init)
    index_t0 = btime.index(t0)
    out = open("me_v"+init+".txt", "w")
    print("center & ME of speed (m/s) & RMSE of speed (m/s) & ME of direction (degree) & RMSE of direction (degree) \\\ ", file=out)
    print("\\hline \\hline", file=out)
    for orig in origin:
        track = np.loadtxt(orig+"/velocity"+init+".txt")
        time = []
        for t in track[:, (0, 1, 2, 3)].astype(np.int32):
            time.append(datetime(*t))
        imax = time.index(tmax)
        print(imax)
        me_v = 0.0
        me_d = 0.0
        rms_v = 0.0
        rms_d = 0.0
        i = 0
        for t in time:
            i = time.index(t)
            if i > imax : break
            ib = btime.index(t)
            v = track[i, 4] - v_bst[ib]
            d = track[i, 5] - d_bst[ib]
            me_v += v
            me_d += d
            rms_v += v**2
            rms_d += d**2
        me_v = me_v / (imax+1)
        me_d = me_d / (imax+1)
        rms_v = np.sqrt(rms_v / (imax+1))    
        rms_d = np.sqrt(rms_d / (imax+1))    
        print("{} & {:7.4f} & {:7.4f} & {:7.4f} & {:7.4f} \\\ ".format(center[orig],\
                 me_v, rms_v, me_d, rms_d), file=out)
    t0 += dt0