import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import sys
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['axes.labelsize']  = 18
plt.rcParams['axes.titlesize']  = 18

def uv2pol(u,v):
    r = np.sqrt(u**2+v**2)
    theta = np.arctan(v / u)
    if theta.ndim > 0:
        for i in range(theta.size):
            if theta[i] < 0.0:
                theta[i] += np.pi
    else:
        if theta < 0.0:
            theta += np.pi
    return r, theta
# main
plot_hod = False
init = sys.argv[1] #yyyymmddhh
iyr  = int(init[0:4])
imo  = int(init[4:6])
idy  = int(init[6:8])
ihr  = int(init[8:])
base  = datetime( iyr,imo,idy,ihr)
#idate = datetime( iyr,imo,idy,ihr)
idate = datetime(2019, 10, 11,12)
vdate = datetime(2019, 10, 12, 0)
tdelta = timedelta(hours=3)
ndate = int((vdate - idate) / tdelta) + 1
print(ndate)
dates = [idate + tdelta * i for i in range(ndate)]
if idate == vdate:
    title = idate.strftime("%Y%m%d%H")
else:
    title = idate.strftime("%Y%m%d%H") + "-" + vdate.strftime("%Y%m%d%H")
t0 = int((idate - base) / tdelta)
offset_era = int((base - datetime(2019, 10, 9, 0)) / tdelta)
t = t0
bstfile = "../pytrack/bst_velocity.txt"
nbst = 0
#sstlist = ["clim", "est", "mgd"]
#exp= {"clim":"clim", "est":"est", "mgd":"mgdsst"}
#suffixes= {"clim":"", "est":"_est", "mgd":"_mgdsst"}
#colors  = {"clim":"red", "est":"blue", "mgd":"tab:green"}
prtblist = ["cntl","en-","en+", "tc"] #, "prtbf"]
exp = {"cntl":"cntl", "en-":"p", "tc":"p2", "en+":"pn", "prtbf":"pf"}
suffixes = {"cntl":"","en-":"+p","tc":"+p2","en+":"+pn","prtbf":"+pf"}
colors = {"cntl":"red","en-":"blue","tc":"tab:green","en+":"cyan","prtbf":"magenta"}
#fig = plt.figure(figsize=(8,8))
#ax = fig.add_subplot(projection="polar")
for date in dates:
    #fig = plt.figure(figsize=(8,8))
    #ax = fig.add_subplot(projection="polar")
    cdate = date.strftime("%Y%m%d%H") 
    print(cdate)
    # best track (only velocity)
    vdata = np.loadtxt(bstfile)
    for i in range(vdata.shape[0]):
        bdate = datetime(int(vdata[i,0]), int(vdata[i,1]), int(vdata[i,2]), int(vdata[i,3]))
        if bdate == date:
            print(f"{i} bdate={bdate.strftime('%Y%m%d%H')}")
            break
    if i < vdata.shape[0]-1:
        vwe = vdata[i,8]
        vns = vdata[i,9]
        if nbst == 0:
            bvwe = 0.0
            bvns = 0.0
        bvwe += vwe
        bvns += vns
        nbst += 1
        #r_v, theta_v = uv2pol(vwe, vns)
        #ax.scatter(theta_v,r_v,s=50,marker='x',color='k',label='best')
    # ERA5
    mfile = f"steer_era5_2019100900/hod_{cdate}.txt"
    sfile = f"steer_era5_2019100900/steer_mean.txt"
    vfile = "../pytrack/velocity_era5.txt"
    data = np.loadtxt(mfile)
    #print(data.shape)
    ndata = data.shape[1]
    level_r = data[0,:]
    u = data[1,:]
    v = data[2,:]
    if t==t0:
        u_r = np.zeros_like(u)
        v_r = np.zeros_like(v)
    u_r += u
    v_r += v
    #r, theta = uv2pol(u,v)
    #ax.plot(theta, r, lw=1, ms=2, marker='o', color='tab:orange')#, label="reanl")
    #for k in range(0,level_r.size,4):
    #    ax.annotate(f"{int(level_r[k])}",xy=(theta[k],r[k]),xycoords='data')
    sdata = np.loadtxt(sfile)
    us = sdata[0,t+offset_era]
    vs = sdata[1,t+offset_era]
    if t==t0:
        us_r = 0.0
        vs_r = 0.0
    us_r += us
    vs_r += vs
    #r_s, theta_s = uv2pol(us, vs)
    #ax.scatter(theta_s,r_s,s=50,marker='^',edgecolors='tab:orange',facecolors='none')
    ##,label="layer mean")
    vdata = np.loadtxt(vfile)
    print(f"ERA5 {int(vdata[t+offset_era,0])} {int(vdata[t+offset_era,1])}"+\
        f" {int(vdata[t+offset_era,2])} {int(vdata[t+offset_era,3])}")
    vwe = vdata[t+offset_era,8]
    vns = vdata[t+offset_era,9]
    if t==t0:
        vwe_r = 0.0
        vns_r = 0.0
    vwe_r += vwe
    vns_r += vns
    #r_v, theta_v = uv2pol(vwe, vns)
    #ax.scatter(theta_v,r_v,s=50,marker='x',color='tab:orange'
    ##,label="TC motion")
    #,label="reanl")
    #ax.tick_params()
    #angle = np.deg2rad(130)
    #ax.legend(loc='lower right',
    #bbox_to_anchor=(.5+np.cos(angle)/2, .5+np.sin(angle)/2))
    #ax.set_title(cdate)
    #fig.savefig(f"hod_{cdate}_era5.png")
    #plt.close()
    if t==t0:
        #nlists = len(sstlist)
        nlists = len(prtblist)
        u_e = np.zeros((nlists,16))
        v_e = np.zeros((nlists,16))
        us_e = np.zeros(nlists)
        vs_e = np.zeros(nlists)
        vwe_e = np.zeros(nlists)
        vns_e = np.zeros(nlists)
    k=0
    #for sst in sstlist:
    for prtb in prtblist:
        #fig = plt.figure(figsize=(8,8))
        #ax = fig.add_subplot(projection="polar")
    #    mfile = f"steer_{exp[sst]}_{init}/hod_{cdate}.txt"
    #    sfile = f"steer_{exp[sst]}_{init}/steer_mean.txt"
    #    vfile = f"../pytrack/velocity{init}_gsm_tl959{suffixes[sst]}.txt"
        mfile = f"steer_{exp[prtb]}_{init}/hod_{cdate}.txt"
        sfile = f"steer_{exp[prtb]}_{init}/steer_mean.txt"
        vfile = f"../pytrack/velocity{init}_gsm_tl479_est{suffixes[prtb]}.txt"
        data = np.loadtxt(mfile)
   #print(data.shape)
        ndata = data.shape[1]
        level_e = data[0,:]
        u = data[1,:]
        v = data[2,:]
        u_e[k,:] += u
        v_e[k,:] += v
        #r, theta = uv2pol(u,v)
        #ax.plot(theta, r, lw=1, ms=2, marker='o', color=colors[sst])#, label=sst)
        #ax.plot(theta, r, lw=1, ms=2, marker='o', color=colors[prtb])#, label=sst)
        #for l in range(0,level_e.size,2):
        #    ax.annotate(f"{int(level_e[l])}",xy=(theta[l],r[l]),
        #    xycoords='data',horizontalalignment='left')
        sdata = np.loadtxt(sfile)
        us = sdata[0,t]
        vs = sdata[1,t]
        us_e[k] += us
        vs_e[k] += vs
        #r_s, theta_s = uv2pol(us, vs)
        #ax.scatter(theta_s,r_s,s=50,marker='^',edgecolors=colors[sst],facecolors='none')
        #ax.scatter(theta_s,r_s,s=50,marker='^',edgecolors=colors[prtb],facecolors='none')
        #,label="layer mean")
        vdata = np.loadtxt(vfile)
        vwe = vdata[t,8]
        vns = vdata[t,9]
        vwe_e[k] += vwe
        vns_e[k] += vns
        #r_v, theta_v = uv2pol(vwe, vns)
        #ax.scatter(theta_v,r_v,s=50,marker='x',color=colors[sst]
        #,label=sst)
        #ax.scatter(theta_v,r_v,s=50,marker='x',color=colors[prtb]
        #,label=prtb)
        #,label="TC motion")
        #ax.tick_params()
        #angle = np.deg2rad(130)
        #ax.legend(loc='lower right',
        #bbox_to_anchor=(.5+np.cos(angle)/2, .5+np.sin(angle)/2))
        #ax.set_title(cdate)
        #fig.savefig(f"hod_{cdate}_{prtb}.png")
        #plt.close()
        k+=1
    #ax.tick_params()
    #angle = np.deg2rad(165)
    #ax.legend(loc='lower right',
    #bbox_to_anchor=(.5+np.cos(angle)/3, .5+np.sin(angle)/3))
    #ax.set_title(cdate)
    #fig.savefig(f"hod_tl959/hod_{cdate}.png")
    #plt.close()
    t += 1
#exit()
level_list = [850, 500, 300]
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(projection="polar")
# best track
bvwe /= nbst
bvns /= nbst
r_v, theta_v = uv2pol(bvwe, bvns)
ax.scatter(theta_v,r_v,s=100,marker='x',color='k',label='Best')
# ERA5
u_r /= ndate
v_r /= ndate
r, theta = uv2pol(u_r,v_r)
if plot_hod:
    fig_h = plt.figure(figsize=(8,8))
    ax_h = fig_h.add_subplot(projection="polar")
    ax_h.plot(theta, r, lw=1, ms=2, marker='o', color='tab:orange')#, label=sst)
    for k in range(level_r.size):
        if int(level_r[k]) % 100 == 0:
            ax_h.annotate(f"{int(level_r[k])}",xy=(theta[k],r[k]),
                xycoords='data',horizontalalignment='left',
                alpha=0.7,zorder=0)
us_r /= ndate
vs_r /= ndate
r_s, theta_s = uv2pol(us_r, vs_r)
ax.scatter(theta_s,r_s,s=100,marker='^',edgecolors='tab:orange',facecolors='none')
if plot_hod:
    ax_h.scatter(theta_s,r_s,s=100,marker='^',edgecolors='tab:orange',facecolors='none'
        ,label='layer mean')
vwe_r /= ndate
vns_r /= ndate
r_v, theta_v = uv2pol(vwe_r, vns_r)
ax.scatter(theta_v,r_v,s=100,marker='x',color='tab:orange',label="Reanl")
if plot_hod:
    ax_h.scatter(theta_v,r_v,s=100,marker='x',color='tab:orange'
        ,label="TC motion")
    ax_h.tick_params()
    #ax.set_ylim(0.0,22.5)
    angle = np.deg2rad(120)
    ax_h.legend(loc='lower right',
    bbox_to_anchor=(.45+np.cos(angle)/2, .45+np.sin(angle)/2))
    ax_h.set_title(title)
    fig_h.savefig(f"hod_{title}_era5.png")
    #plt.close()
# Experiments
k=0
#for sst in sstlist:
for prtb in prtblist:
    u_e[k] /= ndate
    v_e[k] /= ndate
    r, theta = uv2pol(u_e[k],v_e[k])
    if plot_hod:
        fig_h = plt.figure(figsize=(8,8))
        ax_h = fig_h.add_subplot(projection="polar")
        #ax_h.plot(theta, r, lw=1, ms=2, marker='o', color=colors[sst])
        ax_h.plot(theta, r, lw=1, ms=2, marker='o', color=colors[prtb])
        for j in range(level_e.size):
            if int(level_e[j]) % 100 == 0:
                ax_h.annotate(f"{int(level_e[j])}",xy=(theta[j],r[j]),
                xycoords='data',horizontalalignment='left',
                alpha=0.7,zorder=0)
    #ax.plot(theta, r, lw=1, ms=2, marker='o', color=colors[prtb],label=prtb.upper())
    #for j in range(level_e.size):
    #    if int(level_e[j]) in level_list:
    #        ax.annotate(f"{int(level_e[j])}",xy=(theta[j],r[j]),
    #            xycoords='data',horizontalalignment='left',
    #            alpha=0.7,zorder=0)
    us_e[k] /= ndate
    vs_e[k] /= ndate
    r_s, theta_s = uv2pol(us_e[k], vs_e[k])
    #ax.scatter(theta_s,r_s,s=100,marker='^',edgecolors=colors[sst],facecolors='none'
    ax.scatter(theta_s,r_s,s=100,marker='^',edgecolors=colors[prtb],facecolors='none')
    if plot_hod:
        #ax_h.scatter(theta_s,r_s,s=100,marker='^',edgecolors=colors[sst],facecolors='none'
        ax_h.scatter(theta_s,r_s,s=100,marker='^',edgecolors=colors[prtb],facecolors='none'
        ,label='layer mean')
    vwe_e[k] /= ndate
    vns_e[k] /= ndate
    r_v, theta_v = uv2pol(vwe_e[k], vns_e[k])
    #ax.scatter(theta_v,r_v,s=100,marker='x',color=colors[sst],label=sst)
    ax.scatter(theta_v,r_v,s=100,marker='x',color=colors[prtb],label=prtb.upper())
    if plot_hod:
        #ax_h.scatter(theta_v,r_v,s=100,marker='x',color=colors[sst]
        ax_h.scatter(theta_v,r_v,s=100,marker='x',color=colors[prtb]
            ,label='TC motion')
        ax_h.tick_params()
        #ax.set_ylim(0.0,22.5)
        angle = np.deg2rad(120)
        ax_h.legend(loc='lower right',
        bbox_to_anchor=(.45+np.cos(angle)/2, .45+np.sin(angle)/2))
        ax_h.set_title(title)
        #fig_h.savefig(f"hod_tl959/hod_{title}_{sst}.png")
        fig_h.savefig(f"hod_tl479/hod_{title}_{prtb}.png")
        #plt.close()
    k+=1
ax.tick_params()
ax.set_ylim(0.0,8.0)
angle = np.deg2rad(165)
ax.legend(loc='lower right',
    bbox_to_anchor=(.5+np.cos(angle)/3, .5+np.sin(angle)/3))
ax.set_title(title)
fig.savefig(f"hod_{title}.png")