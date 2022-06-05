import numpy as np
import matplotlib.pyplot as plt
import sys
orig = "ukmo"
normalize = True
multi_center = False
if len(sys.argv) > 1:
    orig = sys.argv[1]
orig_list = ["ecmwf","jma","ncep","ukmo"]
member_dict = {"ecmwf":50,"jma":26,"ncep":20,"ukmo":17}
colors      = {"ecmwf":'blue',"jma":'red',"ncep":'green',"ukmo":'yellowgreen'}
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15,15))
finfo_list = []
slp_list = []
for orig in orig_list:
#if True:
    member = member_dict[orig]
    header = f"prtb_{orig}/finfo_0912_{orig}_"
    finfo = np.zeros((2,3,member))
    #if orig == "jma":
    #finfo = np.zeros((2,3,member+1))
    for m in range(member):
        fili = header + f"{m+1:02d}"
        finfo[:,:,m] = np.loadtxt(fili, max_rows=2)
        if finfo[0,1,m] > 180.0:
            finfo[0,1,m] -= 360.0
        if finfo[1,1,m] > 180.0:
            finfo[1,1,m] -= 360.0
    #if orig == "jma":
    #fili = header + "EnSVSA"
    #finfo[:,:,member] = np.loadtxt(fili, max_rows=2)
    #if finfo[0,1,member] > 180.0:
    #    finfo[0,1,member] -= 360.0
    #if finfo[1,1,member] > 180.0:
    #    finfo[1,1,member] -= 360.0
    #print(finfo)
    finfo_list.append(finfo)
    slp = np.loadtxt(f"prtb_{orig}/ridge_{orig}_0912.txt")
    print(slp)
    slp_list.append(slp[:member])
if len(finfo_list) > 1:
    multi_center = True
    color_list = []
    member_list = []
    mem_all = 0
    for i in range(len(finfo_list)):
        orig = orig_list[i]
        member = member_dict[orig]
        color_list += [colors[orig]]*member
        member_list += [m+1 for m in range(member)]
        mem_all += member
    member = mem_all
    finfo = np.concatenate(finfo_list,axis=2)
    print(finfo.shape)
    print(member_list)
    print(color_list)
if len(slp_list) > 1:
    slp = np.concatenate(slp_list)
#fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15,15))
y = slp 
if normalize:
    y = y - y.mean()
    y = y / y.std()
x = finfo[0,0,:] # amplitude of WN0
if normalize:
    x = x - x.mean()
    x = x / x.std()
for m in range(member):
    if multi_center:
        if member_list[m] < 10:
            size=8
        else:
            size=12
        axes[0,0].plot(x[m],y[m],color=color_list[m],\
            marker=r'${:2d}$'.format(member_list[m]),markersize=size)
    else:
        if m < 9:
            size=8
        else:
            size=12
        axes[0,0].plot(x[m],y[m],color=colors[orig],\
            marker=r'${:2d}$'.format(m+1),markersize=size)
#if orig == "jma":
        axes[0,0].plot(x[member],y[member],marker='^',
        color=colors[orig],markersize=15)
xlim = max(abs(x.min()),x.max())
ylim = max(abs(y.min()),y.max())
axes[0,0].set_xlim(-xlim-0.1*x.std(),xlim+0.1*x.std())
axes[0,0].set_ylim(-ylim-0.1*y.std(),ylim+0.1*y.std())
x = finfo[1,1,:] # phase of WN1
if normalize:
    x = x - x.mean()
    x = x / x.std()
for m in range(member):
    if multi_center:
        if member_list[m] < 10:
            size=8
        else:
            size=12
        axes[0,1].plot(x[m],y[m],color=color_list[m],\
            marker=r'${:2d}$'.format(member_list[m]),markersize=size)
    else:
        if m < 9:
            size=8
        else:
            size=12
        axes[0,1].plot(x[m],y[m],color=colors[orig],\
            marker=r'${:2d}$'.format(m+1),markersize=size)
#if orig == "jma":
        axes[0,1].plot(x[member],y[member],marker='^',
        color=colors[orig],markersize=15)
#    axes[0,1].plot(x[member],y[member],marker='^',color='tab:green')
xlim = max(abs(x.min()),x.max())
ylim = max(abs(y.min()),y.max())
axes[0,1].set_xlim(-xlim-0.1*x.std(),xlim+0.1*x.std())
axes[0,1].set_ylim(-ylim-0.1*y.std(),ylim+0.1*y.std())
x = finfo[0,2,:] / (finfo[0,2,:] + finfo[1,2,:])# Variance ratio WN0/(WN0+WN1)
if normalize:
    x = x - x.mean()
    x = x / x.std()
for m in range(member):
    if multi_center:
        if member_list[m] < 10:
            size=8
        else:
            size=12
        axes[1,0].plot(x[m],y[m],color=color_list[m],\
            marker=r'${:2d}$'.format(member_list[m]),markersize=size)
    else:
        if m < 9:
            size=8
        else:
            size=12
        axes[1,0].plot(x[m],y[m],color=colors[orig],\
            marker=r'${:2d}$'.format(m+1),markersize=size)
#if orig == "jma":
        axes[1,0].plot(x[member],y[member],marker='^',
        color=colors[orig],markersize=15)
#    axes[1,0].plot(x[member],y[member],marker='^',color='tab:green')
xlim = max(abs(x.min()),x.max())
ylim = max(abs(y.min()),y.max())
if normalize:
    axes[1,0].set_xlim(-xlim-0.1*x.std(),xlim+0.1*x.std())
else:
    axes[1,0].set_xlim(0.0,1.0)
axes[1,0].set_ylim(-ylim-0.1*y.std(),ylim+0.1*y.std())
#x = finfo[1,2,:] # variance ratio of WN1
#if normalize:
#    x = x - x.mean()
#    x = x / x.std()
#for m in range(member):
#    if m < 9:
#        size=8
#    else:
#        size=12
#    if multi_center:
#        axes[1,1].plot(x[m],y[m],color=color_list[m],\
#            marker=r'${:2d}$'.format(member_list[m]),markersize=size)
#    else:
#        axes[1,1].plot(x[m],y[m],color=colors[orig],\
#            marker=r'${:2d}$'.format(m+1),markersize=size)
#if orig == "jma":
#    axes[1,1].plot(x[member],y[member],marker='^',color='tab:green')
#xlim = max(abs(x.min()),x.max())
#ylim = max(abs(y.min()),y.max())
#axes[1,1].set_xlim(-xlim-0.1*x.std(),xlim+0.1*x.std())
#axes[1,1].set_ylim(-ylim-0.1*y.std(),ylim+0.1*y.std())
axes[0,0].set_xlabel("WN0")
axes[0,1].set_xlabel("Phase of WN1")
axes[1,0].set_xlabel("Variance ratio WN0/(WN0+WN1)")
if not normalize:
    axes[1,0].vlines([0.5],0,1,color='gray',transform=axes[1,0].get_xaxis_transform(),zorder=0)
#axes[1,1].set_xlabel("Variance ratio of WN1")
for ax in axes.flatten():
    ax.hlines(0,1,[0.0],color='gray',transform=ax.get_yaxis_transform(),zorder=0)
    ax.vlines([0.0],0,1,color='gray',transform=ax.get_xaxis_transform(),zorder=0)
    ax.set_ylabel("Perturbation at Ridge")
axes[1,1].remove()
fig.tight_layout()
if multi_center:
    if normalize:
        fig.savefig(f"finfo_scatter_0912_norm.png", dpi=150)
    else:
        fig.savefig(f"finfo_scatter_0912.png", dpi=150)
else:
    if normalize:
        fig.savefig(f"prtb_{orig}/finfo_scatter_0912_{orig}_norm.png", dpi=150)
    else:
        fig.savefig(f"prtb_{orig}/finfo_scatter_0912_{orig}.png", dpi=150)