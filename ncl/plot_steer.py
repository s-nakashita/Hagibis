import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sin, cos, deg2rad
from datetime import datetime, timedelta
import sys

#nlon = 1440 #TL959
nlon = 720 #TL479
lon = np.arange(nlon) * 360.0 / nlon
print(lon)

init = "2019100900"
exp = "clim"
if len(sys.argv) > 1:
    init = sys.argv[1]
if len(sys.argv) > 2:
    exp  = sys.argv[2]
iyr  = int(init[0:4])
imo  = int(init[4:6])
idy  = int(init[6:8])
ihr  = int(init[8:])
base = datetime(iyr, imo, idy, ihr)
end  = datetime(2019, 10, 12, 12)
tdelta = timedelta(hours=3)
date = base
datadir = f"steer_{exp}_{init}/"
mfile = "steer_mean.txt"
dates = ""
wsmean = []
while date <= end:
    cdate = date.strftime("%Y%m%d%H")
    print(cdate)
    hfile = datadir + "hod_" + cdate + ".txt"
    dates = dates + cdate + " "
    #unp = np.loadtxt(datadir + "us_" + cdate + ".txt")
    #print(f"unp min:{np.min(unp)} max:{np.max(unp)}")
    #vnp = np.loadtxt(datadir + "vs_" + cdate + ".txt")
    #print(f"vnp min:{np.min(vnp)} max:{np.max(vnp)}")
    try:
        data = np.loadtxt(datadir+"ushod_"+cdate+".mtx", skiprows=1)
    except FileNotFoundError:
        break
    unp = data[:,1:]
    print(f"unp shape:{unp.shape} min:{np.min(unp)} max:{np.max(unp)}")
    data = np.loadtxt(datadir+"vshod_"+cdate+".mtx", skiprows=1)
    vnp = data[:,1:]
    print(f"vnp shape:{vnp.shape} min:{np.min(vnp)} max:{np.max(vnp)}")
    nlev = unp.shape[1]
    level = np.loadtxt(datadir+"ushod_"+cdate+".mtx", usecols=range(1,nlev+1), max_rows=1)
    print(level)
    i_btm = np.argmin(np.abs(level-850.0))
    i_top = np.argmin(np.abs(level-300.0))
    wgt = np.zeros(i_top - i_btm+1)
    dpall = 850.0 - 300.0
    for i in range(wgt.size-1):
        dp = level[i+i_btm] - level[i+i_btm+1]
        wgt[i] += dp / dpall * 0.5
        wgt[i+1] += dp / dpall * 0.5
    print(wgt)
    print(wgt.sum())
    date += tdelta
    us = np.zeros_like(unp)
    vs = np.zeros_like(vnp)
    for i in range(nlon):
        uvnp = np.array([unp[i],vnp[i]])
        #print(uvnp.shape)
        Arot = np.array(
            [[cos(deg2rad(lon[i])), -sin(deg2rad(lon[i]))],
            [sin(deg2rad(lon[i])), cos(deg2rad(lon[i]))]]
        )
        #print(Arot)
        uvs = np.dot(Arot, uvnp)
        us[i] = uvs[0]
        vs[i] = uvs[1]
    print(f"us shape:{us.shape} min:{np.min(us)} max:{np.max(us)}")
    print(f"vs shape:{vs.shape} min:{np.min(vs)} max:{np.max(vs)}")
    ws = np.sqrt(us**2 + vs**2)
    sdata = np.array([level,us.mean(axis=0),vs.mean(axis=0),ws.mean(axis=0),
    us.std(axis=0),vs.std(axis=0),ws.std(axis=0)])
    np.savetxt(hfile, sdata, fmt="%10.5f")
    us = np.sum(us[:,i_btm:i_top+1] * wgt[None,:], axis=1)
    vs = np.sum(vs[:,i_btm:i_top+1] * wgt[None,:], axis=1)
    ws = np.sum(ws[:,i_btm:i_top+1] * wgt[None,:], axis=1)
    wsmean.append(np.array([us.mean(), vs.mean(), ws.mean(), 
                            us.std(), vs.std(), ws.std()]))
#
#    fig, axs = plt.subplots(2,1,sharex=True,figsize=(12,8))
#    axs[0].plot(lon, unp, color='tab:blue', 
#    marker='x', linewidth=0.0, label='u')
#    axs[0].hlines(unp.mean(), 0, 1, 
#    transform=axs[0].get_yaxis_transform(),color='tab:blue')
#    axs[0].plot(lon, vnp, color='tab:orange',
#    marker='x', linewidth=0.0, label='v')
#    axs[0].hlines(vnp.mean(), 0, 1, 
#    transform=axs[0].get_yaxis_transform(),color='tab:orange')
#    axs[0].set_title('TC coordinates')
#    axs[1].plot(lon,  us, color='tab:blue',
#    marker='x', linewidth=0.0, label='u')
#    axs[1].hlines(us.mean(), 0, 1, 
#    transform=axs[1].get_yaxis_transform(),color='tab:blue')
#    axs[1].plot(lon,  vs, color='tab:orange',
#    marker='x', linewidth=0.0, label='v')
#    axs[1].hlines(vs.mean(), 0, 1, 
#    transform=axs[1].get_yaxis_transform(),color='tab:orange')
#    axs[1].set_title('lonlat coordinates')
#    axs[1].set_xlabel('longitude')
#    for ax in axs:
#        ax.hlines([0.0], 0, 1, transform=ax.get_yaxis_transform(),linestyle='dotted',color='k')
#        ax.set_xticks(lon[0::240])
#        ax.legend()
#    fig.savefig(datadir + "steer_" + cdate + ".png")
#    plt.close()
np.savetxt(datadir + mfile, np.array(wsmean).transpose(), header = dates)