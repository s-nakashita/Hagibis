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
while t0 <= t0max:
    init = t0.strftime("%Y%m%d%H") # -> 2019100600
    year = t0.strftime("%Y")
    print(init)
    index_t0 = btime.index(t0)
    for orig in origin:
        track = np.loadtxt(orig+"/velocity"+init+".txt")
        err = np.loadtxt(orig+"/error"+init+".txt")
        time = []
        for t in track[:, (0, 1, 2, 3)].astype(np.int32):
            time.append(datetime(*t))
        out = open(orig+"/err+v"+init+".txt", "w")
        print("forecast date & error (km) & speed diff. (m/s) & direction diff. (degree)  \\\ ", file=out)
        print("\\hline \\hline", file=out)
        for t in time:
            i = time.index(t)
            ib = btime.index(t)
            ft = err[i, 0]
            e = err[i, 1]
            v = track[i, 4] - v_bst[ib]
            d = track[i, 5] - d_bst[ib]
            print("{} & {:6.3f} & {:6.3f} & {:6.3f} \\\ ".format(t.strftime("%d/%HZ"),\
                 e, v, d), file=out)
    t0 += dt0