import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units as munits
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 16
import sys

cp=1005.7*munits('J/kg/K')
Rd=287.04*munits('J/kg/K')
Lh=2.5104e6*munits('J/kg')
Tr=270.0*munits('K')
pr=1000.0*munits('hPa')
def calc_te(data,base,t,data2=None,levb=925.0,levt=300.0):
    plev = data['level'].sel(level=slice(levb,levt)).values
    plevhPa = plev * munits("hPa")
    #print(plev)
    lon = data['lon'].values
    lat = data['lat'].values
    #print(lat)
    dplev = plev[:-1] - plev[1:]
    wgtlev = np.zeros_like(plev)
    wgtlev[0] = 0.5*dplev[0]/pr.magnitude
    wgtlev[1:-1] = 0.5*dplev[1:]/pr.magnitude + 0.5*dplev[:-1]/pr.magnitude
    wgtlev[-1] = 0.5*dplev[-1]/pr.magnitude
    print(wgtlev,wgtlev.sum())
    wgtlat = np.cos(np.deg2rad(lat))
    area = np.sum(wgtlat)*lon.size

    U = data['U'].sel(time=t,level=slice(levb,levt))
    U = U.values * munits('m/s')
    Ub = base['U'].sel(time=t,level=slice(levb,levt))
    Ub = Ub.values * munits('m/s')
    V = data['V'].sel(time=t,level=slice(levb,levt))
    V = V.values * munits('m/s')
    Vb = base['V'].sel(time=t,level=slice(levb,levt))
    Vb = Vb.values * munits('m/s')
    T = data['T'].sel(time=t,level=slice(levb,levt))
    T = T.values * munits('K')
    Tb = base['T'].sel(time=t,level=slice(levb,levt))
    Tb = Tb.values * munits('K')
    Ps = data['P'].sel(time=t)
    Ps = Ps.values * 1e-2 * munits('hPa')
    Psb = base['P'].sel(time=t)
    Psb = Psb.values * 1e-2 * munits('hPa')
    # calculate Q
    if 'Q' in data:
        Q = data['Q'].sel(time=t,level=slice(levb,levt))
        Q = Q.values * munits('kg/kg')
        Qb = base['Q'].sel(time=t,level=slice(levb,levt))
        Qb = Qb.values * munits('kg/kg')
    else:
        if 'TTD' in data:
            TTD = data['TTD'].sel(time=t,level=slice(levb,levt))
            TTD = TTD.values * munits('degC')
            TTDb = base['TTD'].sel(time=t,level=slice(levb,levt))
            TTDb = TTDb.values * munits('degC')
        elif 'RH' in data:
            RH = data['RH'].sel(time=t,level=slice(levb,levt))
            RH = RH.values * munits('%')
            RHb = base['RH'].sel(time=t,level=slice(levb,levt))
            RHb = RHb.values * munits('%')
            TTD = mpcalc.dewpoint_from_relative_humidity(T,RH)
            TTDb = mpcalc.dewpoint_from_relative_humidity(Tb,RHb)
        Q = mpcalc.specific_humidity_from_dewpoint(plevhPa[:,None,None],TTD)
        Qb = mpcalc.specific_humidity_from_dewpoint(plevhPa[:,None,None],TTDb)
    #print(TTD)
    #print(Q)

    # subtracting base
    U1 = U - Ub
    V1 = V - Vb
    T1 = T - Tb
    Q1 = Q - Qb
    Ps1 = Ps - Psb
    if data2 is not None:
        U = data2['U'].sel(time=t,level=slice(levb,levt))
        U = U.values * munits('m/s')
        V = data2['V'].sel(time=t,level=slice(levb,levt))
        V = V.values * munits('m/s')
        T = data2['T'].sel(time=t,level=slice(levb,levt))
        T = T.values * munits('K')
        Ps = data2['P'].sel(time=t)
        Ps = Ps.values * 1e-2 * munits('hPa')
        # calculate Q
        if 'Q' in data2:
            Q = data2['Q'].sel(time=t,level=slice(levb,levt))
            Q = Q.values * munits('kg/kg')
        else:
            if 'TTD' in data2:
                TTD = data2['TTD'].sel(time=t,level=slice(levb,levt))
                TTD = TTD.values * munits('degC')
            elif 'RH' in data2:
                RH = data2['RH'].sel(time=t,level=slice(levb,levt))
                RH = RH.values * munits('%')
                TTD = mpcalc.dewpoint_from_relative_humidity(T,RH)
            Q = mpcalc.specific_humidity_from_dewpoint(plevhPa[:,None,None],TTD)
        #
        U2 = U - Ub
        V2 = V - Vb
        T2 = T - Tb
        Q2 = Q - Qb
        Ps2 = Ps - Psb
        print(U1.max(),U1.min(),U2.max(),U2.min())
        print(V1.max(),V1.min(),V2.max(),V2.min())
        print(T1.max(),T1.min(),T2.max(),T2.min())
        print(Q1.max(),Q1.min(),Q2.max(),Q2.min())
        print(Ps1.max(),Ps1.min(),Ps2.max(),Ps2.min())
    else:
        U2 = U1
        V2 = V1
        T2 = T1
        Q2 = Q1
        Ps2 = Ps1
        print(U1.max(),U1.min())
        print(V1.max(),V1.min())
        print(T1.max(),T1.min())
        print(Q1.max(),Q1.min())
        print(Ps1.max(),Ps1.min())
    
    ke3d = 0.5*(U1*U2 + V1*V2)
    pe3d = 0.5*(cp*T1*T2/Tr)
    peps = 0.5*(Rd*Tr*(Ps1/pr)*(Ps2/pr))
    lh3d = 0.5*(Lh**2*Q1*Q2/cp/Tr)

    ke = np.sum(ke3d.magnitude * wgtlev[:,None,None] * wgtlat[None,:,None])/area
    pe = (np.sum(pe3d.magnitude * wgtlev[:,None,None] * wgtlat[None,:,None])\
        +np.sum(peps.magnitude * wgtlat[:,None]))/area
    lh = np.sum(lh3d.magnitude * wgtlev[:,None,None] * wgtlat[None,:,None])/area

    mte = ke + pe + lh
    dte = ke + pe

    endict = {'mte':mte,'dte':dte,'ke':ke,'pe':pe,'lh':lh}

    return endict

init = datetime(2019,10,9,12)
cntldir = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est')
prtbdir1 = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p')
prtbdir2 = Path('/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+pn')

figdir = Path(f'GSMTl479L100/{init.strftime("%Y%m%d%H")}')
if not figdir.exists():
    figdir.mkdir(parents=True)

header = 'fcst_p_asia'
cntl = xr.open_dataset(cntldir/f'{header}_{init.strftime("%Y%m%d%H")}.nc')
plus = xr.open_dataset(prtbdir1/f'{header}_{init.strftime("%Y%m%d%H")}.nc')
minus = xr.open_dataset(prtbdir2/f'{header}_{init.strftime("%Y%m%d%H")}.nc')
print(cntl)

time = cntl['time']
en1dict = dict()
en2dict = dict()
magdict = dict()
angdict = dict()
i=0
for t in time:
    print(t)
    endict1 = calc_te(plus,cntl,t)
    endict2 = calc_te(minus,cntl,t)
    endict12 = calc_te(plus,cntl,t,data2=minus)
    if i==0:
        for key in endict1.keys():
            en1dict[key] = []
            en2dict[key] = []
            magdict[key] = []
            angdict[key] = []
    for key in endict1.keys():
        e1 = endict1[key]
        e2 = endict2[key]
        en1dict[key].append(e1)
        en2dict[key].append(e2)
        e12 = endict12[key]
        magratio = np.sqrt(e1 / e2)
        angdiff = np.rad2deg(np.arccos(e12/np.sqrt(e1)/np.sqrt(e2)))
        magdict[key].append(magratio)
        angdict[key].append(angdiff)
    i+=1

for key in magdict.keys():
    en1 = en1dict[key]
    en2 = en2dict[key]
    mag = magdict[key]
    ang = angdict[key]
    fig, axs = plt.subplots(nrows=3,figsize=[8,10],sharex=True,constrained_layout=True)
    axs[0].plot(time,en1,label=r'RIDGE$+$')
    axs[0].plot(time,en2,label=r'RIDGE$-$')
    axs[0].set_ylabel('Perturbation energy\n'+r'[J kg$^{-1}$ m$^{-2}$]')
    axs[1].plot(time,mag)
    axs[2].plot(time,ang)
    axs[1].hlines([1.0],0,1,colors='r',transform=axs[1].get_yaxis_transform())
    axs[2].hlines([180.0],0,1,colors='r',transform=axs[2].get_yaxis_transform())
    axs[1].set_ylabel('magnitude ratio')
    axs[2].set_ylabel('angular difference\n[degree]')
    axs[2].set_xlabel('forecast date')
    plt.setp(axs[2].get_xticklabels(), rotation=30, ha="right")
    fig.suptitle(key.upper())
    axs[0].legend()
    fig.savefig(figdir/f'linearity_{header}_{key}.png',dpi=300)
    plt.show()