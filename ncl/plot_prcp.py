import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import sys
plt.rcParams['legend.fontsize'] = 24
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['axes.labelsize']  = 20
plt.rcParams['axes.titlesize']  = 20
init = sys.argv[1] #yyyymmddhh
iyr  = int(init[0:4])
imo  = int(init[4:6])
idy  = int(init[6:8])
ihr  = int(init[8:])
region = sys.argv[2] #region1 or KN or ST

idate = datetime(2019, 10, 11, 1)
vdate = datetime(2019, 10, 12, 13)
base  = datetime(2019, 10, 11, 0)
tdelta = timedelta(hours=1)

fig, ax1 = plt.subplots(figsize=(8,8))
sstlist = ["clim", "est", "mgd"]
suffixes= {"clim":"", "est":"_est", "mgd":"_mgdsst"}
colors  = {"clim":"red", "est":"blue", "mgd":"tab:green"}
width = 0.25*tdelta
#data
pfile = f"prcp_{region}_{init}.txt"
pdata = np.loadtxt(pfile,skiprows=1)
print(pdata.shape)
time = pdata[:,0] #seconds from 2019-10-11 00:00:00
#print(time)
dates = [base + timedelta(seconds=time[i]) for i in range(time.size)]
#print(dates)
#bar plot
prcp = pdata[:,1]
xaxis = np.array(dates) - 1.5*width
ax1.bar(xaxis, prcp, width=width, 
color='k', alpha=0.5, zorder=1)
for k in range(len(sstlist)):
    sst = sstlist[k]
    prcp = pdata[:,2+k]
    prcp[prcp > 1e20] = np.nan
    xaxis = xaxis + width
    ax1.bar(xaxis, prcp, width=width, 
    color=colors[sst], alpha=0.5, zorder=1)
ax1.set_xticks(dates[5::6])
xlabels = [dates[i].strftime("%m-%d %H") for i in range(len(dates))]
ax1.set_xticklabels(xlabels[5::6])
ax1.set_xticks(dates, minor=True)
ax1.set_xlim(idate, vdate)
ax1.grid(True, which='major', axis='x', zorder=0)
ax1.set_ylabel('Hourly precipitation (mm/hr)')
ax1.set_ylim(0.0, 30.0)
for label in ax1.get_xticklabels():
    label.set_rotation(30)
    label.set_horizontalalignment('right')
#curve plot
ax2 = ax1.twinx()
aprcp = pdata[:,5]
ax2.plot(dates, aprcp, linewidth=3.0, color='k', label='Radar')
for k in range(len(sstlist)):
    sst = sstlist[k]
    aprcp = pdata[:,6+k]
    aprcp[aprcp > 1e20] = np.nan
    ax2.plot(dates, aprcp, linewidth=3.0, color=colors[sst], label=sst.upper())
ax2.set_xlim(idate, vdate)
ax2.set_ylabel('Accumulated precipitation (mm)')
if region == "region1":
    ax2.set_ylim(0.0, 200.0)
elif region == "KN":
    ax2.set_ylim(0.0, 300.0)
elif region == "ST":
    ax2.set_ylim(0.0, 200.0)
ax2.legend(loc='upper left')
fig.tight_layout()
fig.savefig("prcp_"+region+"_"+init+".pdf")
plt.show()
