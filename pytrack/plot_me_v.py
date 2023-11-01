import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 18

t0    = datetime(2019, 10,  7,  0)
tmax  = datetime(2019, 10, 12,  0)
tv    = datetime(2019, 10, 12, 12)
dt    = timedelta(hours=12)
ndate = int((tmax - t0)/dt) + 1
#dates = [t0 + dt*i for i in range(ndate)]
dates = [int((tv - (t0 + dt*i))/timedelta(hours=1)) for i in range(ndate)]
print(dates)

center = ["ecmwf", "jma", "ncep", "ukmo"]
#colors = ["blue" , "red", "green","lightgreen"]
colors = ['dimgray','darkgray','lightgray','white']
values = ["me_v", "rms_v", "me_d", "rms_d"]
titles = {"me_v":"ME speed",
          "rms_v":"RMSE speed",
          "me_d":"ME direction",
          "rms_d":"RMSE direction"}
units = {"me_v":"m/s",
        "rms_v":"m/s",
        "me_d":"degree",
        "rms_d":"degree"}
#fig, axs = plt.subplots(2,2,sharex=False,figsize=(16,10))
stat = np.zeros((4,len(center),len(dates)))
me_v = np.zeros((len(center), len(dates)))
rms_v = np.zeros((len(center), len(dates)))
me_d = np.zeros((len(center), len(dates)))
rms_d = np.zeros((len(center), len(dates)))
xdates = []
i = 0
while t0 <= tmax:
    init = t0.strftime("%Y%m%d%H")
    data = np.loadtxt(f"me_v{init}_vonly.txt")
    stat[0,:,i]  = data[:,0]
    stat[1,:,i] = data[:,1]
    stat[2,:,i]  = data[:,2]
    stat[3,:,i] = data[:,3]
    #me_v[:,i]  = data[:,0]
    #rms_v[:,i] = data[:,1]
    #me_d[:,i]  = data[:,2]
    #rms_d[:,i] = data[:,3]
    xdates.append(t0.strftime("%d-%H"))
    t0 += dt
    i += 1
print(xdates)
for i in range(4):
    val = values[i]
    fig, ax = plt.subplots(figsize=(11,6))
    #width = 0.25*timedelta(hours=6)
    width = 2
    base = np.array(dates) - 1.5*width
    for j in range(len(center)):
        ax.bar(base, stat[i,j], width=width, bottom=0.0, \
            color=colors[j],edgecolor='black', \
            label=center[j].upper())
    #axs[0,0].bar(base, me_v[j], width=width, bottom=0.0, color=colors[j], label=center[j])
    #axs[1,0].bar(base, rms_v[j], width=width, bottom=0.0, color=colors[j], label=center[j])
    #axs[0,1].bar(base, me_d[j], width=width, bottom=0.0, color=colors[j], label=center[j])
    #axs[1,1].bar(base, rms_d[j], width=width, bottom=0.0, color=colors[j], label=center[j])
        base += width
    if val[:2] == "me":
#        ax.hlines([0], dates[0]-dt, dates[-1]+dt, linestyle="solid", color="gray")
#        ax.set_xlim(dates[0]-dt, dates[-1]+dt)
        ax.hlines([0], 0, 1, transform=ax.get_yaxis_transform(), linestyle="solid", color="gray")
#axs[0,0].hlines([0], dates[0], dates[-1], linestyle="solid", color="gray")
#axs[0,1].hlines([0], dates[0], dates[-1], linestyle="solid", color="gray")
#for ax in axs.flatten():
    ax.grid(axis='x')
    ax.legend()
#for ax in axs[1,:]:
    ax.set_xticks(np.arange(132,0,-12))
    ax.set_xlim(138,6)
    ax.set_xlabel('hours from '+tv.strftime("%Y-%m-%d %H00 UTC"))
    #ax.set_xticks(dates)
    #ax.set_xticklabels(xdates)
    #for label in ax.get_xticklabels():
    #    label.set_rotation(30)
    #    label.set_horizontalalignment('right')
    ax.set_title(titles[val],fontsize=24)
    ax.set_ylabel(units[val])
    fig.savefig(f"{val}_tigge_mono.png",dpi=300)
#axs[0,0].set_title("ME speed")
#axs[0,0].set_ylabel("m/s")
#axs[1,0].set_title("RMSE speed")
#axs[1,0].set_ylabel("m/s")
#axs[0,1].set_title("ME direction")
#axs[0,1].set_ylabel("degree")
#axs[1,1].set_title("RMSE direction")
#axs[1,1].set_ylabel("degree")
#fig.savefig("me_v_tigge.png")