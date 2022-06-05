import numpy as np
import matplotlib.pyplot as plt
import sys
orig = "jma"
if len(sys.argv) > 1:
    orig = sys.argv[1]
member_dict = {"ecmwf":50,"jma":26,"ncep":20,"ukmo":17}
member = member_dict[orig]
header = f"prtb_{orig}/finfo_0912_{orig}_"
#if orig == "jma":
finfo = np.zeros((2,3,member+1))
#else:
#    finfo = np.zeros((2,3,member))
for m in range(member):
    fili = header + f"{m+1:02d}"
    finfo[:,:,m] = np.loadtxt(fili, max_rows=2)
    if finfo[0,1,m] > 180.0:
        finfo[0,1,m] -= 360.0
    if finfo[1,1,m] > 180.0:
        finfo[1,1,m] -= 360.0
#if orig == "jma":
fili = header + "EnSVSA"
finfo[:,:,member] = np.loadtxt(fili, max_rows=2)
if finfo[0,1,member] > 180.0:
    finfo[0,1,member] -= 360.0
if finfo[1,1,member] > 180.0:
    finfo[1,1,member] -= 360.0
finfo[:,2,:] *= 100 #%
print(finfo)

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15,10))
#if orig == "jma":
xaxis = np.arange(member+1)+1
#else:
#    xaxis = np.arange(member)+1
width = 0.35
axes[0,1].bar(xaxis,finfo[1,1,:],width=width,color='tab:orange',label='WN1')
axes[0,1].hlines(0,member+2,[0.0],color='gray')
slp = np.loadtxt(f"prtb_{orig}/ridge_{orig}_0912.txt")
axes[1,1].bar(xaxis,slp,width=width,color='tab:green',label='Ridge')
axes[1,1].hlines(0,member+2,[0.0],color='gray')
xaxis2 = xaxis - width*0.5
for j in range(2):
    axes[0,0].bar(xaxis2,finfo[j,0,:],width=width,label=f'WN{j}')
    axes[1,0].bar(xaxis2,finfo[j,2,:],width=width,label=f'WN{j}')
    xaxis2 += width
axes[0,0].hlines(0,member+2,[0.0],color='gray')
axes[0,0].set_title("Amplitude (hPa)")
axes[0,1].set_title("Phase (degree)")
axes[1,0].set_title("Variance (%)")
axes[1,1].set_title("Perturbation (hPa)")
for ax in axes.flatten():
    ax.set_xlabel("member")
    xticks = [f"{i:02d}" for i in range(1,member+1)]
    #if orig == "jma":
    ax.set_xlim(0,member+2)
    xticks += ["mode1"]
    #else:
    #    ax.set_xlim(0,member+1)
    ax.set_xticks(xaxis,labels=xticks)#,minor=True)
    ax.legend()
    for label in ax.get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('right')
fig.tight_layout()
fig.savefig(f"prtb_{orig}/finfo_0912_{orig}.png", dpi=150)