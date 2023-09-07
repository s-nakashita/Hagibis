import numpy as np
from numpy import pi, sin, cos
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import sys
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['axes.labelsize']  = 18
plt.rcParams['axes.titlesize'] = 18
deg2rad = pi / 180.0
init = sys.argv[1]
iyr = int(init[0:4])
imo = int(init[4:6])
idy = int(init[6:8])
ihr = int(init[8:])

plotreanl = False
fig, axall  = plt.subplots(figsize=(12,5)) #All
figwe, axwe = plt.subplots(figsize=(12,5)) #West-East
figns, axns = plt.subplots(figsize=(12,5)) #North-South
# best track
data = np.loadtxt("bst_velocity.txt")
bstlat = []
bstv = []
bstvwe = []
bstvns = []
i = 0
for j in range(data.shape[0]):
    lat = data[j,5]
    v   = data[j,6]
    vwe = data[j,8]
    vns = data[j,9]
    bstlat.append(lat)
    bstv.append(v)
    bstvwe.append(vwe)
    bstvns.append(vns)
axall.plot(bstlat, bstv,  linewidth=3.0, color="black", label="best track")
axwe .plot(bstlat, bstvwe, linewidth=3.0, color="black", label="best track")
axns .plot(bstlat, bstvns, linewidth=3.0, color="black", label="best track")
minlat = np.min(np.array(bstlat))
maxlat = np.max(np.array(bstlat))
## ERA5
if plotreanl:
    data = np.loadtxt("velocity_era5.txt")
    elat = []
    ev   = []
    evwe = []
    evns = []
    i = 0
    k = 0
    for j in range(data.shape[0]):
        lat = data[j,5]
        v   = data[j,6]
        vwe = data[j,8]
        vns = data[j,9]
        elat.append(lat)
        ev.append(v)
        evwe.append(vwe)
        evns.append(vns)
    axall.plot(elat, ev, linewidth=1.0, marker='o', color="orange", label="reanl")
    axwe .plot(elat, evwe, linewidth=1.0, marker='o', color="orange", label="reanl")
    axns .plot(elat, evns, linewidth=1.0, marker='o', color="orange", label="reanl")
    minlat = min(minlat, np.min(np.array(elat)))
    maxlat = max(maxlat, np.max(np.array(elat)))
# experiment
sstlist = ["", "_est", "_mgdsst"]
sstlegend = ["clim", "est", "mgd"]
colors = ["red", "blue", "tab:green"]
#prtblist = ["", "+p", "+p2", "+pn", "+pf"]
#prtblegend = ["cntl", "prtb1", "prtb2", "prtb1_n", "prtbf"]
#colors = ["red", "blue", "tab:green", "cyan", "darkorange"]
minlat = 90.0
maxlat = -90.0
minv   = 999.
minvwe = 999.
minvns = 999.
maxv   = 0.0
maxvwe = 0.0
maxvns = 0.0
for k in range(len(sstlist)):
#for k in range(len(prtblist)):
    vfile = f"velocity{init}_gsm_tl959{sstlist[k]}.txt"
#    vfile = f"velocity{init}_gsm_tl479_est{prtblist[k]}.txt"
    data = np.loadtxt(vfile)
    elat = []
    ev   = []
    evwe = []
    evns = []
    for j in range(data.shape[0]):
        lat = data[j,5]
        v   = data[j,6]
        vwe = data[j,8]
        vns = data[j,9]
        elat.append(lat)
        ev.append  (v)
        evwe.append(vwe)
        evns.append(vns)
    axall.plot(elat,  ev, linewidth=1.0, marker='o', 
    color=colors[k], label=sstlegend[k])
    #color=colors[k], label=prtblegend[k])
    axwe.plot(elat, evwe, linewidth=1.0, marker='o', 
    color=colors[k], label=sstlegend[k])
    #color=colors[k], label=prtblegend[k])
    axns.plot(elat, evns, linewidth=1.0, marker='o', 
    color=colors[k], label=sstlegend[k])
    #color=colors[k], label=prtblegend[k])
    minlat = min(minlat, np.min(np.array(elat)))
    minv   = min(minv,   np.min(np.array(ev  )))
    minvwe = min(minvwe, np.min(np.array(evwe)))
    minvns = min(minvns, np.min(np.array(evns)))
    maxlat = max(maxlat, np.max(np.array(elat)))
    maxv   = max(maxv,   np.max(np.array(ev))  )
    maxvwe = max(maxvwe, np.max(np.array(evwe)))
    maxvns = max(maxvns, np.max(np.array(evns)))
axall.set_ylim(minv-0.5, maxv+0.5)
axwe .set_ylim(minvwe-0.5, maxvwe+0.5)
axns .set_ylim(minvns-0.5, maxvns+0.5)
# (only W-E) plot zero line
start = minlat-0.5
end = maxlat+0.5
axwe.hlines([0], start, end, linestyle="dotted", color="gray")
# rotate labels
for ax in [axall, axwe, axns]:
    ax.legend(ncol=2)
    ax.set_xlabel("latitude")
    ax.set_xlim(start, end)
    #for label in ax.get_xticklabels():
    #    label.set_rotation(30)
    #    label.set_horizontalalignment('right')
axall.set_title("speed (m/s)")
axwe.set_title("W-E speed (m/s)")
axns.set_title("N-S speed (m/s)")
fig.tight_layout()
figwe.tight_layout()
figns.tight_layout()
fig.savefig(f"velocity-lat{init}_gsmtl959.png")
figwe.savefig(f"velocity-lat{init}_gsmtl959_WE.png")
figns.savefig(f"velocity-lat{init}_gsmtl959_NS.png")
plt.close()