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

idate = datetime(2019, 10, 9, 0)
vdate = datetime(2019, 10, 12, 12)
base  = datetime(iyr, imo, idy, ihr)
#tdelta = timedelta(hours=6)

fig, ax = plt.subplots(2,1,figsize=(12,10))
sstlist = ["clim", "est", "mgd"]
suffixes= {"clim":"", "est":"_est", "mgd":"_mgdsst"}
colors  = {"clim":"red", "est":"blue", "mgd":"tab:green"}
for k in range(len(sstlist)):
    sst = sstlist[k]
    flfile= "fllh_avr_core_"+init+suffixes[sst]+".txt"
    fsfile= "flsh_avr_core_"+init+suffixes[sst]+".txt"
    wfile = "rwspd_avr_core_"+init+suffixes[sst]+".txt"
    fldata = np.loadtxt(flfile)
    fldata[0] = np.nan
    print(fldata.shape)
    fsdata = np.loadtxt(fsfile)
    fsdata[0] = np.nan
    print(fsdata.shape)
    wdata = np.loadtxt(wfile)
    print(wdata.shape)
    ndate = fldata.size
    tdelta = timedelta(hours=6)
    width = 0.25*tdelta
    fdates = [base + width*(k-1) + tdelta*i - 0.5*tdelta for i in range(ndate)]
    ax[0].bar(fdates, fldata, width=width, 
    color=colors[sst], edgecolor='white', label=sst.upper(), zorder=1)
    #ax[0].bar(fdates, fsdata, width=width, bottom=fldata,
    #color=colors[sst], edgecolor='white',alpha=0.5, zorder=1)
    ndate = wdata.size
    tdelta = timedelta(hours=3)
    wdates = [base + tdelta*i for i in range(ndate)]
    ax[1].plot(wdates, wdata, linewidth=3.0, color=colors[sst]\
            ,label=sst.upper())
tdelta = timedelta(hours=6)
dates = [base + tdelta*i for i in range(len(fdates))]
ax[0].set_xticks(dates, minor=True)
ax[0].set_xlim(idate, vdate)
ax[0].grid(True, which='both', axis='x', zorder=0)
ax[1].set_xticks(dates, minor=True)
ax[1].set_xlim(idate, vdate)
ax[1].grid(True, which='both', axis='x', zorder=0)
ax[0].legend()
ax[1].legend()
for label in ax[0].get_xticklabels():
    label.set_rotation(30)
    label.set_horizontalalignment('right')
for label in ax[1].get_xticklabels():
    label.set_rotation(30)
    label.set_horizontalalignment('right')
ax[0].set_ylabel(r'Latent heat flux (W/m$^2$)')
#ax[0].set_ylabel(r'Surface heat flux (W/m$^2$)')
ax[1].set_ylabel('10 m radial wind (m/s)')
#ax[0].set_title('average within 200 km radius')
fig.tight_layout()
fig.savefig("fllh+rwspd_avr_core_"+init+".pdf")
plt.show()
