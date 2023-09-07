import numpy as np
from numpy import pi, sin, cos
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import sys
plt.rcParams['legend.fontsize'] = 24
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 20
deg2rad = pi / 180.0
init = sys.argv[1]
iyr = int(init[0:4])
imo = int(init[4:6])
idy = int(init[6:8])
ihr = int(init[8:])

#base = datetime(iyr, imo, idy, ihr)
base = datetime(2019, 10, 9, 0)
vdate = datetime(2019, 10, 12, 12)
n6h = int((datetime(2019, 10, 11, 12) - base)/timedelta(hours=6)) + 1
print(n6h)
bstdates = [base + timedelta(hours=(6*i)) for i in range(n6h)]
base = datetime(2019, 10, 11, 15)
n3h = int((vdate-base)/timedelta(hours=3)) + 1
bstdates = bstdates + [base + timedelta(hours=(3*i)) for i in range(n3h)]
Ndate = len(bstdates)
print(bstdates)

plotreanl = False
fig, axall  = plt.subplots(2,1,sharex=True,figsize=(12,8)) #All
figwe, axwe = plt.subplots(2,1,sharex=True,figsize=(12,8)) #West-East
figns, axns = plt.subplots(2,1,sharex=True,figsize=(12,8)) #North-South
# best track
data = np.loadtxt("bst_velocity.txt")
bstv = []
bstvwe = []
bstvns = []
i = 0
for j in range(data.shape[0]):
    idate = datetime(int(data[j,0]), int(data[j,1]), int(data[j,2]), int(data[j,3]))
    if idate == bstdates[i]:
        v   = data[j,6]
        vwe = data[j,8]
        vns = data[j,9]
        bstv.append(v)
        bstvwe.append(vwe)
        bstvns.append(vns)
        i += 1
    if i == Ndate:
        break
label="BEST"
axall[0].plot(bstdates, bstv,  linewidth=3.0, color="black", label=label)
axwe[0].plot(bstdates, bstvwe, linewidth=3.0, color="black", label=label)
axns[0].plot(bstdates, bstvns, linewidth=3.0, color="black", label=label)
base = datetime(iyr, imo, idy, ihr)
Ndate = int((vdate-base)/timedelta(hours=3))+1
dates = [base + timedelta(hours=(3*i)) for i in range(Ndate)]
print(dates)
## ERA5
if plotreanl:
    width = 0.25*timedelta(hours=3)
    base = np.array(bstdates) - 2*width#*timedelta(hours=6)
    data = np.loadtxt("velocity_era5.txt")
    ev   = []
    evwe = []
    evns = []
    err   = np.zeros(len(bstdates))
    errwe = np.zeros(len(bstdates))
    errns = np.zeros(len(bstdates))
    i = 0
    k = 0
    for j in range(data.shape[0]):
        idate = datetime(int(data[j,0]), int(data[j,1]), int(data[j,2]), int(data[j,3]))
        if idate == dates[i]:
            v   = data[j,6]
            vwe = data[j,8]
            vns = data[j,9]
            ev.append(v)
            evwe.append(vwe)
            evns.append(vns)
            i += 1
        if idate == bstdates[k]:
            err[k:]   += v   - bstv[k]
            errwe[k:] += vwe - bstvwe[k]
            errns[k:] += vns - bstvns[k]
            k += 1
        if i == Ndate:
            break
    axall[0].plot(dates, ev, linewidth=3.0, marker='o', color="orange", label="reanl")
    axwe [0].plot(dates, evwe, linewidth=3.0, marker='o', color="orange", label="reanl")
    axns [0].plot(dates, evns, linewidth=3.0, marker='o', color="orange", label="reanl")
    axall[1].bar(base, err, width=width, bottom=0.0, color="orange", label="reanl")
    axwe [1].bar(base, errwe, width=width, bottom=0.0, color="orange", label="reanl")
    axns [1].bar(base, errns, width=width, bottom=0.0, color="orange", label="reanl")
    base = base + width#*timedelta(hours=6)
else:
    width = 0.3*timedelta(hours=3)
    base = np.array(bstdates) - 1.0*width#*timedelta(hours=6)

# experiment
sstlist = ["", "_est", "_mgdsst"]
sstlegend = ["clim", "est", "mgd"]
colors = ["red", "blue", "tab:green"]
#prtblist = ["", "+p", "+p2", "+pn", "+pf"]
#prtblegend = ["cntl", "prtb1", "prtb2", "prtb1_n", "prtbf"]
#colors = ["red", "blue", "tab:green", "cyan", "darkorange"]
for k in range(len(sstlist)):
#for k in range(len(prtblist)):
    vfile = f"velocity{init}_gsm_tl959{sstlist[k]}.txt"
#    vfile = f"velocity{init}_gsm_tl479_est{prtblist[k]}.txt"
    data = np.loadtxt(vfile)
    ev   = []
    evwe = []
    evns = []
    err   = np.zeros(len(bstdates))
    errwe = np.zeros(len(bstdates))
    errns = np.zeros(len(bstdates))
    i = 0
    if ihr == 0:
        l = 0
    else:
        l = 2
    for j in range(data.shape[0]):
        idate = datetime(int(data[j,0]), int(data[j,1]), int(data[j,2]), int(data[j,3]))
        if idate == dates[i]:
            v   = data[j,6]
            vwe = data[j,8]
            vns = data[j,9]
            ev.append  (v)
            evwe.append(vwe)
            evns.append(vns)
            i += 1
        if idate == bstdates[l]:
            err[l:]   += v   - bstv[l]
            errwe[l:] += vwe - bstvwe[l]
            errns[l:] += vns - bstvns[l]
            l += 1
        if i == Ndate:
            break
    label=sstlegend[k].upper()
    #label=prtblegend[k].upper()
    axall[0].plot(dates,  ev, linewidth=3.0, #marker='o', 
    color=colors[k], label=label)
    axwe[0].plot(dates, evwe, linewidth=3.0, #marker='o', 
    color=colors[k], label=label)
    axns[0].plot(dates, evns, linewidth=3.0, #marker='o', 
    color=colors[k], label=label)
    axall[1].bar(base,  err, width=width, bottom=0.0, 
    color=colors[k], label=label)
    axwe[1].bar(base, errwe, width=width, bottom=0.0, 
    color=colors[k], label=label)
    axns[1].bar(base, errns, width=width, bottom=0.0, 
    color=colors[k], label=label)
    base = base + width#*timedelta(hours=6)
# (only W-E and bars) plot zero line
dt = timedelta(hours=3)
idate = bstdates[0]
start = bstdates[0]-dt
end = dates[-1]+dt
axall[1].hlines([0], start, end, linestyle="solid", color="gray")
axwe[0].hlines([0], start, end, linestyle="dotted", color="gray")
axwe[1].hlines([0], start, end, linestyle="solid", color="gray")
axns[1].hlines([0], start, end, linestyle="solid", color="gray")
axwe[0].set_ylim(-5,10)
axns[0].set_ylim(2,14)
# rotate labels
for ax in axall:
    ax.set_xlim(idate - dt, vdate + dt)
    ax.set_xticks(dates, minor=True)
    ax.legend(ncol=2)
    ax.grid(True, which='major', axis='x', linestyle='dotted', zorder=0)
    for label in ax.get_xticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment('right')
for ax in axwe:
    ax.set_xlim(idate - dt, vdate + dt)
    ax.set_xticks(dates, minor=True)
    ax.legend(ncol=2)
    ax.grid(True, which='major', axis='x', linestyle='dotted', zorder=0)
    for label in ax.get_xticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment('right')
for ax in axns:
    ax.set_xlim(idate - dt, vdate + dt)
    ax.set_xticks(dates, minor=True)
    ax.legend(ncol=2)
    ax.grid(True, which='major', axis='x', linestyle='dotted', zorder=0)
    for label in ax.get_xticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment('right')
axall[0].set_title("speed (m/s)")
axall[1].set_title("accumulated speed error (m/s)")
axwe[0].set_title("W-E speed (m/s)")
axwe[1].set_title("accumulated W-E speed error (m/s)")
axns[0].set_title("N-S speed (m/s)")
axns[1].set_title("accumulated N-S speed error (m/s)")
fig.savefig(f"velocity{init}_gsmtl959.pdf")
figwe.savefig(f"velocity{init}_gsmtl959_WE.pdf")
figns.savefig(f"velocity{init}_gsmtl959_NS.pdf")
plt.show()
plt.close()
