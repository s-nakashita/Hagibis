import numpy as np
import matplotlib.pyplot as plt
nlev=3
orig='ukmo'
if orig == 'jma':
    nlev = 6
etype='dry'
TE = np.loadtxt(f"TE-{etype}-m1-{orig}-2019100912-2019101212_nlev{nlev}.txt")
KE = np.loadtxt(f"KE-{etype}-m1-{orig}-2019100912-2019101212_nlev{nlev}.txt")
PE = np.loadtxt(f"PE-{etype}-m1-{orig}-2019100912-2019101212_nlev{nlev}.txt")
if etype != 'dry':
    LE = np.loadtxt(f"LE-{etype}-m1-{orig}-2019100912-2019101212_nlev{nlev}.txt")
else:
    LE = np.zeros_like(TE)
print(TE[:,1])
print(KE[:,1]+PE[:,1]+LE[:,1])

time = TE[:,0]
vtime = np.arange(time.size)[time==72.0]
it_v = vtime[0]
print(it_v)
if etype == 'moist':
    if orig == 'ecmwf':
        sigma = 73.9661
    if orig == 'jma':
        sigma = 65.7825
    if orig == 'ncep':
        sigma = 42.1070
elif etype == "wmoist":
    if orig == 'jma':
        sigma = 40.8240
elif etype == 'dry':
    if orig == 'jma':
        sigma = 43.3248
    if orig == 'ukmo':
        sigma = 37.7293
print(f"eigen value ={sigma**2:16.8e}")
maxTE = (sigma**2)*TE[0,1]
evrate = TE[it_v,1] / TE[0,1]
print(f"evolution rate ={evrate:16.8e}")
print(f"ratio ={evrate/(sigma**2):16.8e}")

fig, ax = plt.subplots(figsize=(10,6))
width=5.0
#ax.plot(time, TE[:,1], label='total')
ax.bar(time, KE[:,1], width, label='kinetic')
ax.bar(time, PE[:,1], width, bottom=KE[:,1], label='potential')
if etype == 'moist' or etype == 'wmoist':
    ax.bar(time, LE[:,1], width, bottom=KE[:,1]+PE[:,1], label='latent heat')
#ax.plot([time[0],time[it_v]], [TE[0,1],maxTE], linestyle='dashed', color='gray')
ax.plot(time, maxTE*np.ones(time.size), linestyle='dashed', color='gray')
ax.vlines(time[it_v], 0, 1, transform=ax.get_xaxis_transform(), colors='r', linestyle='dashdot')
ax.set_yscale('log')
ax.set_ylim(1e-2,1e3)
ax.set_xlim(time[0],time[-1])
ax.set_xticks(time[0::2])
ax.set_xticks(time, minor=True)
ax.set_xlabel('forecast time[h]')
ax.set_ylabel('J/kg')
ax.legend()
fig.savefig(f'tevol-{etype}_{orig}_nlev{nlev}.png',bbox_inches='tight',dpi=300)