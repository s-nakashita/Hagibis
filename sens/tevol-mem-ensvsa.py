import numpy as np
import matplotlib.pyplot as plt
nlev=3
orig='ukmo'
if orig == 'jma':
    nlev = 6
etype='dry'
TE = np.loadtxt(f"TE-{etype}-m1-{orig}-2019100912-2019101212_nlev{nlev}.txt")
print(TE.shape)

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
maxTE = TE[0,1]*(sigma**2)

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(time, TE[:,2], color='gray', label='member')
for i in range(3,TE.shape[1]):
    ax.plot(time, TE[:,i], color='gray')
ax.plot(time, TE[:,1], label='EnSVSA')
ax.vlines(time[it_v], 0, 1, transform=ax.get_xaxis_transform(), colors='r', linestyle='dashdot')
ax.set_yscale('log')
ax.set_ylim(1e-2,1e3)
ax.set_xlim(time[0],time[-1])
ax.set_xticks(time[0::2])
ax.set_xticks(time, minor=True)
ax.set_xlabel('forecast time[h]')
ax.set_ylabel('J/kg')
ax.legend()
fig.savefig(f'tevol-{etype}-mem_{orig}_nlev{nlev}.png',bbox_inches='tight',dpi=300)