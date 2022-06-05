import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import sys
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['axes.labelsize']  = 18
plt.rcParams['axes.titlesize']  = 18

def correlation(u,v):
    var1 = np.sum((u - u.mean())**2)
    var2 = np.sum((v - v.mean())**2)
    cov  = np.sum((u-u.mean())*(v-v.mean()))
    corr = cov / np.sqrt(var1 * var2)
    return corr
# main
init = sys.argv[1] #yyyymmdd
iyr  = int(init[0:4])
imo  = int(init[4:6])
idy  = int(init[6:8])
ihr  = int(init[8:])

idate = datetime(2019, 10,  9,  0)
vdate = datetime(2019, 10, 12, 12)
tdelta = timedelta(hours=3)
ndate = int((vdate - idate) / tdelta) + 1
print(ndate)
dates = [idate + tdelta * i for i in range(ndate)]
markdate = vdate + 0.5*tdelta

fig, axs = plt.subplots(2,1,sharex=True,figsize=(12,8))
# ERA5
base = idate
mfile = f"steer_era5_2019100900/steer_mean.txt"
vfile = f"../pytrack/velocity_era5.txt"
data = np.loadtxt(mfile)
print(data.shape)
ndata = data.shape[1]
u = data[0,:]
v = data[1,:]
w = data[2,:]
us = data[3,:]
vs = data[4,:]
ws = data[5,:]
vdata = np.loadtxt(vfile)
vwe = vdata[:ndata,8]
vns = vdata[:ndata,9]
# correlation
corr_u = correlation(u, vwe)
print(f"correlation between W-E steering flow and velocity = {corr_u}")
corr_v = correlation(v, vns)
print(f"correlation between S-N steering flow and velocity = {corr_v}")
#
xaxis = [base + tdelta * i for i in range(data.shape[1])]
axs[0].plot(xaxis, u, color='tab:orange',linewidth=1.0, marker='^',label='reanl')
#axs[0].plot(xaxis, vwe, color='tab:orange',linewidth=1.0, marker='o',label='TC velocity')
axs[1].plot(xaxis, v, color='tab:orange',linewidth=1.0, marker='^',label='reanl')
#axs[1].plot(xaxis, vns, color='tab:orange',linewidth=1.0, marker='o',label='reanl')
#axs[2].plot(xaxis, w, color='tab:orange',linewidth=1.0, marker='^',label='reanl')
#axs[0].errorbar(xaxis, u, yerr=us, color='tab:orange',label='reanl')
#axs[1].errorbar(xaxis, v, yerr=vs, color='tab:orange',label='reanl')
#axs[2].errorbar(xaxis, w, yerr=ws, color='tab:orange',label='reanl')
#sstlist = ["clim", "est", "mgd"]
#exp= {"clim":"clim", "est":"est", "mgd":"mgdsst"}
#suffixes= {"clim":"", "est":"_est", "mgd":"_mgdsst"}
#colors  = {"clim":"red", "est":"blue", "mgd":"tab:green"}
prtblist = ["cntl","prtb1","prtb2", "prtb1_n", "prtbf"]
exp = {"cntl":"cntl", "prtb1":"p", "prtb2":"p2", "prtb1_n":"pn", "prtbf":"pf"}
suffixes = {"cntl":"","prtb1":"+p","prtb2":"+p2","prtb1_n":"+pn","prtbf":"+pf"}
colors = {"cntl":"red","prtb1":"blue","prtb2":"tab:green","prtb1_n":"cyan","prtbf":"magenta"}
#style = {0:"dashed",12:"solid"}
style = {0:"solid",12:"solid"}
#for ihr in [0, 12]:
base  = datetime( iyr,imo,idy,ihr)
#for sst in sstlist:
#    print(sst)
for prtb in prtblist:
    print(prtb)
#        mfile = f"steer_{exp[sst]}_{init}{ihr:02d}/steer_mean.txt"
    #mfile = f"steer_{exp[sst]}_{init}/steer_mean.txt"
    mfile = f"steer_{exp[prtb]}_{init}/steer_mean.txt"
    #vfile = f"../pytrack/velocity{init}_gsm_tl959{suffixes[sst]}.txt"
    vfile = f"../pytrack/velocity{init}_gsm_tl479_est{suffixes[prtb]}.txt"
    data = np.loadtxt(mfile)
    print(data.shape)
    ndata = data.shape[1]
    u = data[0,:]
    v = data[1,:]
    w = data[2,:]
    us = data[3,:]
    vs = data[4,:]
    ws = data[5,:]
    vdata = np.loadtxt(vfile)
    vwe = vdata[:ndata,8]
    vns = vdata[:ndata,9]
    # correlation
    corr_u = correlation(u, vwe)
    print(f"correlation between W-E steering flow and velocity = {corr_u}")
    corr_v = correlation(v, vns)
    print(f"correlation between S-N steering flow and velocity = {corr_v}")
#
    xaxis = [base + tdelta * i for i in range(data.shape[1])]
#        if ihr == 12:
    #axs[0].plot(xaxis, u, color=colors[sst], linewidth=1.0, marker='^',
    #linestyle=style[ihr], label=sst)
    axs[0].plot(xaxis, u, color=colors[prtb], linewidth=1.0, marker='^',
    #axs[0].errorbar(xaxis, u, yerr=us, color=colors[sst], 
    linestyle=style[ihr], label=prtb)
    #axs[0].plot(xaxis, vwe, color=colors[sst], linewidth=1.0, marker='o',
    #linestyle=style[ihr], label="TC velocity")
    #axs[1].plot(xaxis, v, color=colors[sst], linewidth=1.0, marker='^',
    #linestyle=style[ihr], label=sst)
    axs[1].plot(xaxis, v, color=colors[prtb], linewidth=1.0, marker='^',
    #axs[1].errorbar(xaxis, v, yerr=vs, color=colors[sst], 
    linestyle=style[ihr], label=prtb)
    #axs[1].plot(xaxis, vns, color=colors[sst], linewidth=1.0, marker='o',
    #linestyle=style[ihr], label="TC velocity")
    #axs[2].plot(xaxis, w, color=colors[sst], linewidth=1.0, marker='^',
    ##axs[2].errorbar(xaxis, w, yerr=ws, color=colors[sst], 
    #linestyle=style[ihr], label=sst)
#        else:
#            axs[0].plot(xaxis, u, color=colors[sst], linewidth=1.0, marker='^', fillstyle="none",
#            #axs[0].errorbar(xaxis, u, yerr=us, color=colors[sst], 
#            linestyle=style[ihr])
#            axs[1].plot(xaxis, v, color=colors[sst], linewidth=1.0, marker='^', fillstyle="none",
#            #axs[1].errorbar(xaxis, v, yerr=vs, color=colors[sst], 
#            linestyle=style[ihr])
#            axs[2].plot(xaxis, w, color=colors[sst], linewidth=1.0, marker='^', fillstyle="none",
#            #axs[2].errorbar(xaxis, w, yerr=ws, color=colors[sst], 
#            linestyle=style[ihr])
axs[0].set_title('W-E')
#axs[0].set_ylim(-5,10)
axs[1].set_title('S-N')
#axs[1].set_ylim(2,14)
#axs[2].set_title('Speed')
for ax in axs:
    ax.set_xlim(idate - tdelta, vdate + tdelta)
    ax.set_xticks(dates, minor=True)
    ax.grid(True, which='major', axis='x', linestyle='dotted', zorder=0)
for label in axs[1].get_xticklabels():
    label.set_rotation(30)
    label.set_horizontalalignment('right')
axs[0].hlines([0.0], 0, 1, transform=axs[0].get_yaxis_transform(),color='gray',zorder=0)
axs[0].legend()
fig.tight_layout()
fig.savefig("steer_mean_"+init+"_tl479.png")