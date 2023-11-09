import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import sys
plt.rcParams['legend.fontsize'] = 24
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['axes.labelsize']  = 20
plt.rcParams['axes.titlesize']  = 20

def uv2pol(u,v):
    r = np.sqrt(u**2+v**2)
    if u.ndim > 0:
        theta = np.zeros_like(u)
        for i in range(theta.size):
            if abs(u[i]) < 1e-13:
                if v[i] > 0.0:
                    theta[i] = np.pi / 2.
                else:
                    theta[i] = 3.*np.pi / 2.
            else:
                tan = v[i] / u[i]
                theta[i] = np.arctan(tan)
                if theta[i] < 0.0:
                    if v[i] > 0.0:
                        theta[i] += np.pi
                    else:
                        theta[i] += 2*np.pi 
                else:
                    if v[i] < 0.0:
                        theta[i] += np.pi
    else:
        if abs(u) < 1e-13:
            if v > 0.0:
                theta = np.pi / 2.
            else:
                theta = 3.*np.pi / 2.
        else:
            tan = v / u
            theta = np.arctan(tan)
            if theta < 0.0:
                if v > 0.0:
                    theta += np.pi
                else:
                    theta += 2*np.pi 
            else:
                if v < 0.0:
                    theta += np.pi
    return r, theta
# main
plot_hod = True
plot_reanl = False
init = sys.argv[1] #yyyymmddhh
iyr  = int(init[0:4])
imo  = int(init[4:6])
idy  = int(init[6:8])
ihr  = int(init[8:])
base  = datetime( iyr,imo,idy,ihr)
start = init
if len(sys.argv)>2:
    start = sys.argv[2]
iyr  = int(start[0:4])
imo  = int(start[4:6])
idy  = int(start[6:8])
ihr  = int(start[8:])
idate = datetime(iyr,imo,idy,ihr)
vdate = idate + timedelta(hours=24)
#idate = datetime(2019, 10,  9,12)
#vdate = datetime(2019, 10, 10,12)
#idate = datetime(2019, 10, 10, 12)
#vdate = datetime(2019, 10, 11, 12)
#idate = datetime(2019, 10, 11, 12)
#vdate = datetime(2019, 10, 12, 12)
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
#
exptype = "prtb" # sst or prtb
if len(sys.argv) > 3:
    exptype = sys.argv[3]
#
if exptype == "sst":
    tl = 959
    elist = ["est", "clim", "mgd"]
    exp= {"clim":"clim", "est":"est", "mgd":"mgdsst"}
    suffixes= {"clim":"", "est":"_est", "mgd":"_mgdsst"}
    colors  = {"clim":"red", "est":"blue", "mgd":"tab:green"}
elif exptype == "prtb":
    tl = 479
    elist = ["cntl","ridge$+$","ridge$-$"]#,"tc",  "full"]
    exp = {"cntl":"cntl", "ridge$+$":"p", "tc":"p2", "ridge$-$":"pn", "full":"pf"}
    suffixes = {"cntl":"_est","ridge$+$":"_est+p","tc":"_est+p2","ridge$-$":"_est+pn","full":"_est+pf"}
    #colors = {"cntl":"red","ridge$+$":"blue","tc":"tab:green","ridge$-$":"green","full":"magenta"}
    colors = {"cntl":"darkgray","ridge$+$":"darkgray","tc":"darkgray","ridge$-$":"darkgray","full":"darkgray"}
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
    if plot_reanl:
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
        nlists = len(elist)
        u_e = np.zeros((nlists,16))
        v_e = np.zeros((nlists,16))
        us_e = np.zeros(nlists)
        vs_e = np.zeros(nlists)
        vwe_e = np.zeros(nlists)
        vns_e = np.zeros(nlists)
    k=0
    for e in elist:
        #fig = plt.figure(figsize=(8,8))
        #ax = fig.add_subplot(projection="polar")
        mfile = f"steer_{exp[e]}_{init}/hod_{cdate}.txt"
        sfile = f"steer_{exp[e]}_{init}/steer_mean.txt"
        vfile = f"../pytrack/velocity{init}_gsm_tl{tl}{suffixes[e]}.txt"
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
level_list = [1000,850, 500, 300, 200]
mark_list = ['o','s','+']
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(projection="polar")
fig_v = plt.figure(figsize=(8,8))
ax_v = fig_v.add_subplot(projection="polar")
# best track
bvwe /= nbst
bvns /= nbst
r_v, theta_v = uv2pol(bvwe, bvns)
ax_v.scatter(theta_v,r_v,s=100,marker='x',\
        linewidths=2.0,color='k',label='Best')
ax_v.annotate("",
        xytext=(theta_v,r_v),xycoords='data',
        xy=(0,0),textcoords='data',
        horizontalalignment='left',
        arrowprops=dict(arrowstyle="<|-",color='k')
        )
if plot_reanl:
    # ERA5
    u_r /= ndate
    v_r /= ndate
    r, theta = uv2pol(u_r,v_r)
    if plot_hod:
        fig_h = plt.figure(figsize=(8,8))
        ax_h = fig_h.add_subplot(projection="polar")
        ax_h.plot(theta, r, lw=2, ms=2, marker='o', color='tab:orange')#, label=sst)
        for k in range(level_r.size):
            if int(level_r[k]) % 100 == 0:
                ax_h.annotate(f"{int(level_r[k])}",xy=(theta[k],r[k]),
                    xycoords='data',horizontalalignment='left',
                    alpha=0.7,zorder=0)
    us_r /= ndate
    vs_r /= ndate
    r_s, theta_s = uv2pol(us_r, vs_r)
    ax_v.scatter(theta_s,r_s,s=100,marker='^',edgecolors='tab:orange',facecolors='none')
    if plot_hod:
        ax_h.scatter(theta_s,r_s,s=100,marker='^',edgecolors='tab:orange',facecolors='none'
            ,label='layer mean')
    vwe_r /= ndate
    vns_r /= ndate
    r_v, theta_v = uv2pol(vwe_r, vns_r)
    ax_v.scatter(theta_v,r_v,s=100,marker='x',color='tab:orange',label="Reanl")
    if plot_hod:
        ax_h.scatter(theta_v,r_v,s=100,marker='x',color='tab:orange'
        ,label="TC motion")
        ax_h.tick_params()
        ax_h.set_xticks(np.deg2rad(np.array([0.0,90.0,180.0,270.0])))
        ax_h.set_xticklabels(['E','N','W','S'])
        #ax.set_ylim(0.0,22.5)
        angle = np.deg2rad(120)
        ax_h.legend(loc='lower right',
        bbox_to_anchor=(.45+np.cos(angle)/2, .45+np.sin(angle)/2))
        ax_h.set_title(title)
        fig_h.savefig(f"hod_{title}_era5.png")
        plt.close(fig=fig_h)
# Experiments
k=0
for e in elist:
    u_e[k] /= ndate
    v_e[k] /= ndate
    r, theta = uv2pol(u_e[k],v_e[k])
    if plot_hod:
        fig_h = plt.figure(figsize=(8,8))
        ax_h = fig_h.add_subplot(projection="polar")
        ax_h.plot(theta, r, lw=2, ms=2, marker='o', color=colors[e])
        for j in range(level_e.size):
            #if int(level_e[j]) % 100 == 0:
            if int(level_e[j]) in level_list:
                ax_h.annotate(f"{int(level_e[j])}",xy=(theta[j],r[j]),
                xycoords='data',horizontalalignment='left',#alpha=0.7,
                size=16,zorder=0)
    ax.plot(theta, r, lw=2, color=colors[e],label=e.upper())
    imrk=0
    for j in range(level_e.size):
        if int(level_e[j]) in level_list:
            ax.annotate(f"{int(level_e[j])}",xy=(theta[j],r[j]),
                xycoords='data',horizontalalignment='left',
                color=colors[e],size=12,
                zorder=0)
            ax.plot(theta[j],r[j],lw=0,ms=4,marker='o',color=colors[e])
            #imrk+=1
    us_e[k] /= ndate
    vs_e[k] /= ndate
    r_s, theta_s = uv2pol(us_e[k], vs_e[k])
    ### steering flow
    #ax.scatter(theta_s,r_s,s=150,marker='^'\
    #    ,edgecolors=colors[e],facecolors='none',\
    #    linewidths=2.0)
    #if plot_hod:
    #    ax_h.scatter(theta_s,r_s,s=100,marker='^',\
    #        edgecolors=colors[e],facecolors='none',\
    #        linewidths=2.0,label='layer mean')
    vwe_e[k] /= ndate
    vns_e[k] /= ndate
    r_v, theta_v = uv2pol(vwe_e[k], vns_e[k])
    ax.scatter(theta_v,r_v,s=150,marker='x',\
        linewidths=2.0,color=colors[e])
    ax_v.scatter(theta_v,r_v,s=100,marker='x',\
        linewidths=2.0,color=colors[e],label=e.upper())
    ax_v.annotate("",
        xytext=(theta_v,r_v),xycoords='data',
        xy=(0,0),textcoords='data',
        horizontalalignment='left',
        arrowprops=dict(arrowstyle="<|-",color=colors[e])
        )
    if plot_hod:
        deg=90.0-theta_v*180.0/np.pi
        ax_h.scatter(theta_v,r_v,s=100,marker='x',\
            color=colors[e],linewidths=2.0,\
            label='TC motion\n'+f'{r_v:.1f} m/s, {deg:.1f} deg.')
        ax_h.annotate("",
        xytext=(theta_v,r_v),xycoords='data',
        xy=(0,0),textcoords='data',
        horizontalalignment='left',
        arrowprops=dict(arrowstyle="<|-",color='k')
        )
        ax_h.tick_params()
        ax_h.set_xticks(np.deg2rad(np.array([0.0,90.0,180.0,270.0])))
        ax_h.set_xticklabels(['E','N','W','S'])
        if idate.strftime("%Y%m%d%H") == "2019100912":
            ax_h.set_ylim(0.0,8.0)
        elif idate.strftime("%Y%m%d%H") == "2019101012":
            ax_h.set_ylim(0.0,12.0)
        elif idate.strftime("%Y%m%d%H") == "2019101112":
            ax_h.set_ylim(0.0,22.5)
        angle = np.deg2rad(270)
        ax_h.legend(loc='lower center',
        bbox_to_anchor=(np.cos(angle)/2+.5, np.sin(angle)/2+.7))
        ax_h.grid(True)
        ax_h.set_title(f"{e.upper()} {title}")
        fig_h.savefig(f"hod_tl{tl}/hod_mono_{title}_{e}.png",dpi=300)
        fig_h.savefig(f"hod_tl{tl}/hod_mono_{title}_{e}.pdf")
        plt.close(fig=fig_h)
    k+=1
for ax1 in [ax,ax_v]:
    ax1.tick_params()
    ax1.set_xticks(np.deg2rad(np.array([0.0,90.0,180.0,270.0])))
    ax1.set_xticklabels(['E','N','W','S'])
    angle = np.deg2rad(270)
    ax1.legend(loc='lower center',#ncol=2,
        bbox_to_anchor=(np.cos(angle)/2+.5, np.sin(angle)/2+.7))
    ax1.set_title(title)
    ax1.grid(True)
#ax_v.set_ylim(0.0,8.0)
#fig.savefig(f"hod_tl{tl}/hod_mono_{title}.png",dpi=300)
#fig_v.savefig(f"hod_tl{tl}/vel_mono_{title}.png",dpi=300)
#fig.savefig(f"hod_tl{tl}/hod_mono_{title}.pdf")
#fig_v.savefig(f"hod_tl{tl}/vel_mono_{title}.pdf")
#plt.show()