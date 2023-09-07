import numpy as np 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path
from datetime import datetime, timedelta
import sys
import re
plt.rcParams['font.size'] = 18

exps = ['cntl','ridge$+$','ridge$-$']
suffixes = {'cntl':'','ridge$+$':'+p','ridge$-$':'+pn'}
markers = {'cntl':'s','ridge$+$':'o','ridge$-$':'^'}
#colors = {'cntl':'red','ridge$+$':'blue','ridge$-$':'green'}
colors = {'cntl':'dimgray','ridge$+$':'dimgray','ridge$-$':'dimgray'}
styles = {'cntl':'solid','ridge$+$':'dotted','ridge$-$':'dashed'}

sdate = datetime(2019,10,9,12)
edate = datetime(2019,10,9,12)
dt = timedelta(hours=24)
if sdate==edate:
    cdate=sdate.strftime("%Y%m%d%H")
else:
    cdate=sdate.strftime("%Y%m%d%H")+"-"+edate.strftime("%Y%m%d%H")

## best track
lbst = True
yyyy = sdate.strftime("%Y")
fbst = 'bst_hagibis.txt'
try:
    bsttrack = np.loadtxt(fbst)
except FileNotFoundError:
    print(f"not found {fbst}")
    frpd = Path(f'/Volumes/dandelion/data/bsttrack/{yyyy}/rpd{yyyy[2:]}{tcnum:02d}.txt')
    try:
        bsttrack = np.loadtxt(frpd)
    except FileNotFoundError:
        lbst=False
        print(f"not found {frpd}")

# track
fig = plt.figure(figsize=(8,8),constrained_layout=True)
ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
ax.coastlines()
lonw = 130
lone = 148
lats = 15
latn = 45
dlon = dlat = 5
ax.set_xticks(list(range(lonw,lone+dlon,dlon)),crs=ccrs.PlateCarree())
ax.set_yticks(list(range(lats,latn+dlat,dlat)),crs=ccrs.PlateCarree())
ax.set_extent([lonw,lone,lats,latn], ccrs.PlateCarree())
ax.gridlines()

# mslp
fig2, ax2 = plt.subplots(figsize=(10,6),constrained_layout=True)

if lbst:
    lonbst = bsttrack[:,4]
    latbst = bsttrack[:,5]
    slpbst = bsttrack[:,6]
    lonsel = []
    latsel = []
    btime = []
    xticks = []
    texts = []
    for i in range(bsttrack.shape[0]):
        tmp=datetime(int(bsttrack[i,0]),int(bsttrack[i,1]),\
            int(bsttrack[i,2]),int(bsttrack[i,3]))
        btime.append(tmp)
        if int(bsttrack[i,3])==0:
            xticks.append(tmp)
        if tmp==datetime(2019,10,12,12):
            lonsel.append(lonbst[i]);latsel.append(latbst[i])
            texts.append(f"{int(bsttrack[i,2])}")
    ax.plot(lonbst,latbst,c='k',lw=4.0,\
        transform=ccrs.PlateCarree(),label='best track')
    ax.plot(lonsel,latsel,marker='*',lw=0.0,\
        markerfacecolor='w',markeredgecolor='k',markersize=15.0,
        transform=ccrs.PlateCarree())
    #for lon1,lat1,text in zip(lonsel,latsel,texts):
    #    if lon1 < lonw or lon1 > lone or lat1 < lats or lat1 > latn: continue
    #    ax.text(lon1,lat1,text,{'ha':'center','va':'center','c':'k','size':8},\
    #        transform=ccrs.PlateCarree())
    ax2.plot(btime,slpbst,c='k',lw=3.0,label='best track')
    ax2.set_xticks(xticks)

cmap = plt.get_cmap('tab10')
icol = 0
style = 'solid'
suffix = ''
while sdate <= edate:
  for exp in exps:
    #for wdir, style in zip(['rsm2msm9_gfsz'],\
    #    ['solid','dashed']):
    #for wdir, style in zip(['DATA/gfs/rda'],['solid','dashed']):
    tfile = f"track{sdate.strftime('%Y%m%d%H')}_gsm_tl479_est{suffixes[exp]}.txt"
    try:
        track = np.loadtxt(tfile)
    except FileNotFoundError:
        print(f"not found {tfile}")
        sdate+=dt
        continue
    label=exp.upper()
    time = []
    lonsel = []
    latsel = []
    texts = []
    for i in range(track.shape[0]):
        tmp=datetime(int(track[i,0]),int(track[i,1]),\
            int(track[i,2]),int(track[i,3]))
        time.append(tmp)
        if tmp==datetime(2019,10,12,12):
            lonsel.append(track[i,4]);latsel.append(track[i,5])
            texts.append(f"{int(track[i,2])}")
    lonc = track[:,4]
    latc = track[:,5]
    slpc = track[:,6] * 1.0e2 #Pa
    mew=1.5
    ax.plot(lonc,latc,c=colors[exp],lw=2.0,ls=styles[exp],\
        marker=markers[exp],ms=8.0,\
        fillstyle='none',mec=colors[exp],mew=mew,\
        transform=ccrs.PlateCarree())
    mew=3.0
    ax.plot(lonsel,latsel,c=colors[exp],lw=2.0,ls=styles[exp],\
        marker=markers[exp],ms=12.0,\
        fillstyle='none',mec=colors[exp],mew=mew,\
        transform=ccrs.PlateCarree(),label=label)
    ax2.plot(time,slpc,c=colors[exp],lw=3.0,ls=styles[exp],\
        marker=markers[exp],ms=8.0,\
        fillstyle='none',mec=colors[exp],mew=mew,\
        label=label)
    #for lon1,lat1,text in zip(lonsel,latsel,texts):
    #    if lon1 < lonw or lon1 > lone or lat1 < lats or lat1 > latn: continue
    #    ax.text(lon1,lat1,text,{'ha':'center','va':'center','c':'k','size':8},\
    #    transform=ccrs.PlateCarree())
    icol+=1
  sdate+=dt
ax2.grid()
ax.legend()
ax2.legend()
plt.setp(ax2.get_xticklabels(),rotation=30,ha="right")
#fig.savefig(f'track{cdate}_tigge.png',dpi=300)
#fig2.savefig(f'mslp{cdate}_tigge.png',dpi=300)
fig.savefig(f'track{cdate}_gsmtl479_prtb1.png',dpi=300)
fig2.savefig(f'mslp{cdate}_gsmtl479_prtb1.png',dpi=300)
plt.show()
