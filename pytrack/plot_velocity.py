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
    xlabels = [btime[i].strftime("%m%d%H") for i in range(index_t0, index_tmax)]
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    plt.grid(True)
    ax1.plot(btime[index_t0:index_tmax], v_bst[index_t0:index_tmax], color="black", label="best")
    ax2.plot(btime[index_t0:index_tmax], d_bst[index_t0:index_tmax], color="black", label="best")
    for orig in origin:
        track = np.loadtxt(orig+"/velocity"+init+".txt")
        time = []
        for t in track[:, (0, 1, 2, 3)].astype(np.int32):
            time.append(datetime(*t))
        tmax = time.index(datetime(2019, 10, 13, 0))
        v = track[:, 4]
        d = track[:, 5]
        ax1.plot(time[:tmax], v[:tmax], linestyle="None", marker="o", color=mcolor[orig], label=center[orig])
        ax2.plot(time[:tmax], d[:tmax], linestyle="None", marker="o", color=mcolor[orig], label=center[orig])
    ax1.set_title("speed "+init)
    ax2.set_title("direction "+init)
    ax2.invert_yaxis()
    #ax1.set_xticklabels(xlabels,rotation=45, ha='right')
    #ax2.set_xticklabels(xlabels,rotation=45, ha='right')
    fig1.autofmt_xdate(rotation=45)
    fig2.autofmt_xdate(rotation=45)
    ax1.legend()
    ax2.legend()
    fig1.savefig("v"+init+".png")
    fig2.savefig("d"+init+".png")
    t0 += dt0