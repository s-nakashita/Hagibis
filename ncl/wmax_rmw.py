import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import sys
convert_knot = {
    0:0,5:3,10:5,15:8,20:10,25:13,30:15,34:17,
    35:18,40:20,45:23,48:25,50:25,55:30,60:30,64:33,
    65:35,70:35,75:40,80:40,85:45,90:45,95:50,100:50,
    105:55,110:55,115:60,120:60,125:65,130:65,135:70,
    140:70,145:75,150:75,155:80,160:80,165:85,170:85,
    175:90,180:95,185:95,190:100,195:100,200:105,
    205:105,210:110,215:110,220:115,225:115,230:120,
    235:120,235:120,240:125,245:125,250:130,255:130,
    260:135,265:135,270:140,275:140,280:145,285:145,
    290:150,295:150,300:155
}
convert_mile = {
    0:0,0.1:0.2,0.3:0.5,0.5:1,1:2,5:10,10:20,15:30,20:40,25:50,
    30:60,35:60,40:70,45:80,50:90,55:100,60:110,65:120,
    70:130,75:140,80:150,85:160,90:170,95:180,100:190,
    110:200,120:220,130:240,135:250,140:260,150:280,
    160:300,170:310,180:330,190:350,200:370,210:390,
    215:400,220:410,230:430,240:440,250:460,260:480,
    270:500,280:520,290:540,295:550,300:560,325:600,
    350:650,375:700,400:750,425:800,450:850,475:900,
    500:950,550:1000,600:1100,650:1200,700:1300,750:1400,
    800:1500,850:1600,900:1700,950:1800,1000:1900,
    1100:2000,1200:2200,1300:2400,1400:2600,1500:2800,1600:3000
}
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['axes.labelsize']  = 18
plt.rcParams['axes.titlesize']  = 18
fig, axes = \
    plt.subplots(nrows=2,ncols=1,sharex=True,figsize=(11,11))
idate = datetime(2019, 10,  9, 12)
vdate = datetime(2019, 10, 13,  0)
# year mon day hour lon lat slp ws r(ws>25) r(ws>15)
bst = np.loadtxt("../pytrack/bst_hagibis_all.txt")
btime = []
for t in bst[:, (0,1,2,3)].astype(np.int32):
    btime.append(datetime(*t))
wspd = np.zeros(len(btime))
radi = np.zeros(len(btime))
for l in range(len(btime)):
    wspd[l] = convert_knot[bst[l,7].astype(np.int32)]
    radi[l] = convert_mile[bst[l,8].astype(np.int32)]
axes[0].plot(btime, wspd, linewidth=3.0, color='k')
axes[1].plot(btime, radi, linewidth=3.0, color='k')
prtblist = ["cntl", "en-", "en+", "tc"]#,  "prtbf"]
suffixes = {"cntl":"", "en-":"+p", "tc":"+p2", "en+":"+pn", "prtbf":"+pf"}
colors = {"cntl":"red","en-":"blue","tc":"tab:green","en+":"blue","prtbf":"magenta"}
styles = {"cntl":"solid","en-":"solid","tc":"solid","en+":"dashed","prtbf":"solid"}
init = "2019100912"
for prtb in prtblist:
    wfile = "wmax_" + init + "_"+prtb+".txt"
    rfile = "rmw_" + init + "_"+prtb+".txt"
    wdata = np.loadtxt(wfile)
    rdata = np.loadtxt(rfile)
    ndate = wdata.size
    time = []
    tdelta = timedelta(hours=6)
    time = [idate + tdelta*i for i in range(ndate)]
    axes[0].plot(time, wdata, linewidth=2.0, linestyle=styles[prtb], 
        #color=colors[sst], label=sst)
        color=colors[prtb], label=prtb.upper())
    axes[1].plot(time, rdata, linewidth=2.0, linestyle=styles[prtb], 
        #color=colors[sst], label=sst)
        color=colors[prtb], label=prtb.upper())
for ax in axes:
    ax.set_xticks(time[::2])
    ax.set_xticks(time, minor=True)
    xlabels = [time[i].strftime("%m-%d %H") for i in range(len(time))]
    ax.set_xticklabels(xlabels[::2])
    ax.grid(True, which='both', axis='x', zorder=0)
    ax.set_xlim(idate, vdate)
    ax.legend()
    for label in ax.get_xticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment('right')
axes[0].set_ylabel('Max wind speed (m/s)')
axes[0].set_ylim(20.0, 60.0)
axes[1].set_ylabel('Radius of Wspd > 25 m/s (km)')
axes[1].set_ylim(50.0,400.0)
fig.tight_layout()
#fig.savefig("slpc_gsmtl959_0900.png")
fig.savefig("wspd_gsmtl479_0912.png")