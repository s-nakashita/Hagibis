from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['legend.fontsize'] = 24
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['axes.labelsize']  = 20
plt.rcParams['axes.titlesize']  = 20
fig, ax = plt.subplots(figsize=(11,6))
idate = datetime(2019, 10, 9, 12)
vdate = datetime(2019, 10, 13, 0)
# year mon day hour lon   lat  slp
# 2019 10  4   18   164.4 15.7 1008.00
bst = np.loadtxt("bst_hagibis.txt")
btime = []
for t in bst[:, (0, 1, 2, 3)].astype(np.int32):
    btime.append(datetime(*t))
slp = np.zeros(len(btime))
for l in range(len(btime)):
    slp[l] = bst[l, 6]*0.01 # Pa->hPa
ax.plot(btime, slp, linewidth=3.0, color='k')
sstlist = ["clim", "est", "mgd"]
suffixes = {"clim":"", "est":"_est", "mgd":"_mgdsst"}
colors = {"clim":"red", "est":"blue", "mgd":"green"}
prtblist = ["cntl", "ridge$+$", "ridge$-$"]#, "tc",  "full"]
#prtblist = ["cntl", "tc",  "full"]
suffixes = {"cntl":"", "ridge$+$":"+p", "tc":"+p2", "ridge$-$":"+pn", "full":"+pf"}
colors = {"cntl":"red","ridge$+$":"blue","tc":"tab:green","ridge$-$":"green","full":"magenta"}
styles = {"cntl":"solid","ridge$+$":"solid","tc":"solid","ridge$-$":"solid","full":"solid"}
dt0 = timedelta(hours = 12)
#for sst in sstlist:
#    t0 = datetime(2019, 10, 9, 0)
#    t0max = datetime(2019, 10, 9, 0)
for prtb in prtblist:
    t0 = datetime(2019, 10, 9, 12)
    t0max = datetime(2019, 10, 9, 12)
    while t0 <= t0max:
        init = t0.strftime("%Y%m%d%H")
        print(init)
        yyyy = t0.strftime("%Y")
        hh   = t0.strftime("%H")
        #tfile = "track" + init + "_gsm_tl959"+suffixes[sst]+".txt"
        tfile = "track" + init + "_gsm_tl479_est"+suffixes[prtb]+".txt"
        track = np.loadtxt(tfile)
        time = []
        for t in track[:, (0, 1, 2, 3)].astype(np.int32):
            time.append(datetime(*t))
        slp = np.zeros(len(time))
        for l in range(len(time)):
            slp[l] = track[l, 6]
        #if hh == "00":
        #    style = "dashed"
        #    ax.plot(time, slp, linewidth=2.0, linestyle=style, 
        #color=colors[sst])
        #else:
        #style = "solid"
        ax.plot(time, slp, linewidth=3.0, 
        #color=colors[sst], label=sst.upper())
        linestyle=styles[prtb], color=colors[prtb], label=prtb.upper())
        t0 += dt0
bst = np.loadtxt("track_era5.txt")
btime = []
for t in bst[:, (0, 1, 2, 3)].astype(np.int32):
    btime.append(datetime(*t))
slp = np.zeros(len(btime))
for l in range(len(btime)):
    slp[l] = bst[l, 6]*0.01 # Pa->hPa
#ax.plot(btime, slp, linewidth=2.0, color='darkorange', label="Reanl")

#tfile = "track" + init + "_gsm_tl959_est.txt"
tfile = "track" + init + "_gsm_tl479_est.txt"
track = np.loadtxt(tfile)
time = []
for t in track[:, (0, 1, 2, 3)].astype(np.int32):
    time.append(datetime(*t))
ax.set_xticks(time[::4])
ax.set_xticks(time[::2], minor=True)
xlabels = [time[i].strftime("%m-%d %H") for i in range(len(time))]
ax.set_xticklabels(xlabels[::4])
ax.grid(True, which='both', axis='x', zorder=0)
ax.set_xlim(idate, vdate)
ax.legend()
for label in ax.get_xticklabels():
    label.set_rotation(30)
    label.set_horizontalalignment('right')
ax.set_ylabel('central pressure (hPa)')
ax.set_ylim(910.0, 980.0)
fig.tight_layout()
#fig.savefig("slpc_gsmtl959_0900.pdf")
fig.savefig("slpc_gsmtl479_0912_prtb1.pdf")
fig.savefig("slpc_gsmtl479_0912_prtb1.png",dpi=300)
plt.show()
