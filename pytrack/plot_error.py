import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.dates as mdates
plt.rcParams['font.size'] = 18

idate = datetime(2019,10,7,12)
edate = datetime(2019,10,12,12)

origs = ['ecmf','rjtd','kwbc','egrr']
centers = {'ecmf':'ecmwf','rjtd':'jma','kwbc':'ncep','egrr':'ukmo'}
#markers = {'ecmf':'s','rjtd':'o','kwbc':'^','egrr':'p'}
#colors = {'ecmf':'blue','rjtd':'red','kwbc':'green','egrr':'lime'}
colors = {'ecmf':'dimgray','rjtd':'darkgray','kwbc':'lightgray','egrr':'white'}
#hatches = {'ecmf':'--','rjtd':'+','kwbc':'x','egrr':'//'}

fig, ax = plt.subplots(figsize=[8,6],constrained_layout=True)
# 4 centers
width = timedelta(hours=2.5)
distances = {'ecmf':width*(-1.5),'rjtd':width*(-0.5),'kwbc':width*0.5,'egrr':width*1.5}
for orig in origs:
    center = centers[orig]
    err = pd.read_csv(f'error-{orig}.txt',sep=' ', header=None,\
    names=('init','error'))
    xaxis = np.array([datetime.strptime(err.iloc[i,0],'%Y-%m-%dT%H') for i in range(err.shape[0])])
    xaxis += distances[orig]
    print(xaxis)
    rects = ax.bar(xaxis,err['error'],width=width,\
        color=colors[orig],edgecolor='black',#hatch=hatches[orig],
        label=center.upper())
    #ax.bar_label(rects, labels=[center[0].upper() for i in range(xaxis.size)],
    #padding=3, fontsize=16)
    #mew=3.0
    #ax.plot(xaxis,err['error'],\
    #marker=markers[orig],lw=2.0,ms=12.0,\
    #fillstyle='none',mec=colors[orig],mew=mew,\
    #c=colors[orig],label=centers[orig].upper())
# JMA 2019 mean error
jmamean = pd.read_csv('error-mean.txt',sep=' ', header=None,\
    names=('init','error'))
print(jmamean)
xaxis = [datetime.strptime(jmamean.iloc[i,0],'%Y-%m-%dT%H') for i in range(jmamean.shape[0])]
print(xaxis)
ax.plot(xaxis,jmamean['error'],\
    lw=4.0,c='black',label='JMA 2019 mean')
ax.set_ylabel('positional error [km]')
ax.set_xlabel('initial date')
ax.xaxis.set_minor_locator(mdates.HourLocator(interval=12))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H"))
for label in ax.get_xticklabels():
    label.set_rotation(30)
    label.set_horizontalalignment('right')
ax.grid(axis='x',zorder=0)
ax.legend()
fig.savefig('error-tigge.png',dpi=300)
plt.show()