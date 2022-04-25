import numpy as np
import matplotlib.pyplot as plt

latc=35
file=f"_avr_latc{latc}.mtx"
exp="dTE"
out=f"ensvsa-{exp}_prof_latc{latc}.png"
#file="_avr_target.mtx"
#out=f"ensvsa-{exp}_prof_target.png"
#out=f"ensvsa-{exp}_prof.png"
t_range = range(8,13)
data = np.loadtxt("ke"+file, skiprows=1)
ft = data[:,0]
print(ft)
KE = data[:,1:]
print(KE.shape)
nlev = KE.shape[1]-1
plev = np.loadtxt("ke"+file, usecols=range(1,nlev+1), max_rows=1)
print(plev)
p0 = 1000.0
zaxis = -np.log(plev/p0)
fig, ax = plt.subplots(figsize=(8,4))
if exp == "dTE":
    en_list = ["ke","pe"]
else:
    en_list = ["ke","pe","le"]
labels = {"ke":"Kinetic","pe":"Potential","le":"Latent"}
styles = {"ke":"dashed","pe":"dashdot","le":"dotted"}
cmap = plt.get_cmap('tab10')
te = np.zeros((4,ft.size,nlev))
te_int = np.zeros(ft.size)
labels = []
for d in t_range:
    if d < 4:
        latc=15
        file=f"_avr_latc{latc}.mtx"
        labels.append(f"{int(ft[d])}H,lat={latc}")
    elif d < 8:
        latc=25
        file=f"_avr_latc{latc}.mtx"
        labels.append(f"{int(ft[d])}H,lat={latc}")
    #elif d < 12:
    else:
        latc=35
        file=f"_avr_latc{latc}.mtx"
        labels.append(f"{int(ft[d])}H,lat={latc}")
    #else:
    #    file="_avr_target.mtx"
    #    labels.append(f"{int(ft[d])}H,target")
    j=1
    for en in en_list:
        data = np.loadtxt(en+file, \
        usecols=range(1,nlev+1),skiprows=1)
        data_int = np.loadtxt(en+file, \
        usecols=(nlev+1),skiprows=1)
    #i = 0
    #for d in [0,ft.size-1]:
        #if en == "ke":
        #    ax.plot(data[d],zaxis,linestyle=styles[en],color=cmap(i)\
        #    #    ,label=labels[en])
        #        ,label=f"{int(ft[d])}H")
        #else:
        #ax.plot(data[d]/data_int[d],zaxis,linestyle=styles[en],color=cmap(i))
        #i += 1
        te[j,d,:] = data[d,:]
        te[0,d,:] += data[d,:]
        te_int[d] += data_int[d]
        j+=1
i = 0
for d in t_range:
    ax.plot(te[0,d]/te_int[d],zaxis,color=cmap(i),label=labels[i])
    j=1
    for en in en_list:
        ax.plot(te[j,d]/te_int[d],zaxis,color=cmap(i),linestyle=styles[en])
        j+=1
    i += 1
ax.set_ylabel(r'$z^*$')
ax.set_ylim(zaxis[-1],zaxis[0])
#ax.set_xscale("log")
ax.set_xlim(0.0,1.0)
ax.set_xlabel("normalized energy")
ax.legend()
if exp == "TE":
    ax.set_title("solid:Total, dashed:Kinetic, dashdot:Potential, dotted:Latent")
else:
    ax.set_title("solid:Total, dashed:Kinetic, dashdot:Potential")
fig.savefig(out,bbox_inches='tight',dpi=300)