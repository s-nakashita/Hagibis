from librotate import pi, dtheta, mask_lonlat, rotate_lonlat1d
import numpy as np
import matplotlib.pyplot as plt

lonc = 137.5 
latc = 25.0 
dcolat = 12.5 

lon = np.arange(360) 
lat = -90.0 + np.arange(180)
msk = mask_lonlat(dcolat, lonc, latc, lon, lat)
nmsk = np.sum(msk.astype(np.int))
print(nmsk)
mlon = np.zeros(nmsk)
mlat = np.zeros_like(mlon)
print(mlon.size)
k = 0
for j in range(lat.size):
    for i in range(lon.size):
        if msk[j, i]:
            mlon[k] = lon[i]
            mlat[k] = lat[j]
            k += 1
print(mlon.min(),mlon.max())
print(mlat.min(),mlat.max())
lonnp, latnp = rotate_lonlat1d(lonc, latc, mlon, mlat, -1)

fig, ax = plt.subplots(2)
ax[0].scatter(mlon,mlat)
ax[0].scatter(lonc, latc, marker="*")
#ax.set_xlim([lon.min(),lon.max()])
ax[0].set_xlim([90.0,180.0])
#ax.set_ylim([lat.min(),lat.max()])
ax[0].set_ylim([-30.0,60.0])
ax[0].set_aspect("equal")
ax[0].set_xlabel("longitude")
ax[0].set_ylabel("latitude")
ax[1].scatter(lonnp, latnp, s=10)
ax[1].set_xlim([lon.min(),lon.max()])
ax[1].set_ylim([60.0,lat.max()])
ax[1].set_aspect("equal")
ax[1].set_xlabel("longitude")
ax[1].set_ylabel("latitude")
plt.show()
