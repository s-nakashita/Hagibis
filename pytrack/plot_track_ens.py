import numpy as np 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path
from datetime import datetime, timedelta
import sys
import re
plt.rcParams['font.size'] = 18

#origs = ['ecmf','rjtd','kwbc','egrr']
origs = ['rjtd']
centers = {'ecmf':'ecmwf','rjtd':'jma','kwbc':'ncep','egrr':'ukmo'}
markers = {'ecmf':'s','rjtd':'o','kwbc':'^','egrr':'p'}
#colors = {'ecmf':'blue','rjtd':'red','kwbc':'green','egrr':'lime'}
colors = {'ecmf':'black','rjtd':'black','kwbc':'dimgray','egrr':'dimgray'}
member = {'ecmf':50,'rjtd':26,'kwbc':20,'egrr':17}
color_mem = 'dimgray'

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
    ax.plot(lonbst,latbst,c='k',lw=4.0,ls='dashed',\
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
  for orig in origs:
    #for wdir, style in zip(['rsm2msm9_gfsz'],\
    #    ['solid','dashed']):
    #for wdir, style in zip(['DATA/gfs/rda'],['solid','dashed']):
    datadir = Path(f'{orig}')
    tfile = datadir / f"track{sdate.strftime('%Y%m%d%H')}.txt"
    try:
        track = np.loadtxt(tfile)
    except FileNotFoundError:
        print(f"not found {tfile}")
        sdate+=dt
        continue
    label=centers[orig].upper()
    time = []
    lonsel = []
    latsel = []
    texts = []
    for i in range(track.shape[0]):
        tmp=datetime(int(track[i,0]),int(track[i,1]),\
            int(track[i,2]),int(track[i,3]))
        time.append(tmp)
        if tmp==datetime(2019,10,10,12):
            lonsel.append(track[i,4]);latsel.append(track[i,5])
        if tmp==datetime(2019,10,12,12):
            lonsel.append(track[i,4]);latsel.append(track[i,5])
            texts.append(f"{int(track[i,2])}")
            break
    lonc = track[:len(time),4]
    latc = track[:len(time),5]
    slpc = track[:len(time),6]
    mew=1.5
    ax.plot(lonc,latc,c=colors[orig],lw=2.0,ls=style,\
        marker=markers[orig],ms=6.0,\
        fillstyle='none',mec=colors[orig],mew=mew,\
        transform=ccrs.PlateCarree())
    mew=3.0
    ax.plot(lonsel,latsel,c=colors[orig],lw=0.0,ls=style,\
        marker=markers[orig],ms=12.0,\
        fillstyle='none',mec=colors[orig],mew=mew,\
        transform=ccrs.PlateCarree(),label=label)
    ax2.plot(time,slpc,c=colors[orig],lw=3.0,ls=style,\
        marker=markers[orig],ms=8.0,\
        fillstyle='none',mec=colors[orig],mew=mew,\
        label=label)
    #for lon1,lat1,text in zip(lonsel,latsel,texts):
    #    if lon1 < lonw or lon1 > lone or lat1 < lats or lat1 > latn: continue
    #    ax.text(lon1,lat1,text,{'ha':'center','va':'center','c':'k','size':8},\
    #    transform=ccrs.PlateCarree())
    ## member
    mem = member[orig]
    loncmean = np.zeros_like(lonc)
    latcmean = np.zeros_like(latc)
    slpcmean = np.zeros_like(slpc)
    for m in range(1,mem+1):
        datadir = Path(f'{centers[orig]}')
        tfile = datadir / f"track{sdate.strftime('%Y%m%d%H')}_{m:02d}.txt"
        try:
            track = np.loadtxt(tfile)
        except FileNotFoundError:
            print(f"not found {tfile}")
            continue
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
                break
        lonc = track[:len(time),4]
        latc = track[:len(time),5]
        slpc = track[:len(time),6]
        loncmean = loncmean + lonc
        latcmean = latcmean + latc
        slpcmean = slpcmean + slpc
        ax.plot(lonc,latc,c=color_mem,lw=1.0,ls='dotted',\
            transform=ccrs.PlateCarree(),zorder=0)
        mew=1.0
        ax.plot(lonsel,latsel,c=colors[orig],lw=0.0,ls=style,\
            marker=markers[orig],ms=8.0,\
            fillstyle='none',mec=colors[orig],mew=mew,\
            transform=ccrs.PlateCarree(),zorder=0)
        ax2.plot(time,slpc,c=color_mem,lw=3.0,ls=style,zorder=0)
#    loncmean /= float(mem)
#    latcmean /= float(mem)
#    slpcmean /= float(mem)
#    ax.plot(loncmean,latcmean,c='gray',lw=2.0,ls=style,\
#        marker=markers[orig],ms=6.0,\
#        fillstyle='none',mec='gray',mew=mew,\
#        transform=ccrs.PlateCarree())
#    ax2.plot(time,slpc,c='gray',lw=2.0,ls=style,\
#        marker=markers[orig],ms=8.0,\
#        fillstyle='none',mec='gray',mew=mew)
    icol+=1
  sdate+=dt
ax2.grid()
ax.legend()
ax2.legend()
plt.setp(ax2.get_xticklabels(),rotation=30,ha="right")
#figname=f'track{cdate}_ens_tigge'
#figname2=f'mslp{cdate}_ens_tigge'
figname=f'track{cdate}_ens_{centers[orig]}'
figname2=f'mslp{cdate}_ens_{centers[orig]}'
fig.savefig(f'{figname}.png',dpi=300)
fig2.savefig(f'{figname2}.png',dpi=300)
fig.savefig(f'{figname}.pdf')
fig2.savefig(f'{figname2}.pdf')
plt.show()
