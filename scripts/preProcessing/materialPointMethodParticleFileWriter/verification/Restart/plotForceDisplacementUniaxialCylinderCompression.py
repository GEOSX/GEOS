# -*- coding: utf-8 -*-
"""
Created on Mon March 11th 03/11/2024
@author: appleton
Supplemental plot nominal force displacement 
"""
import numpy as np                   # math stuff
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


# parent directory of reaction file path:
runLocation='__PARENT_DIRECTORY_OF_RESULTS__'

files=[
'__SIMULATION_SUBDIRECTORY__/',
]  
labels=[
'__SIMULATION_LABEL__',
]

#initialize plots
evenly_spaced_interval = np.linspace(0, 1, len(files))
colors = [cm.rainbow(x) for x in evenly_spaced_interval]
fig, (ax1, ax2, ax3) = plt.subplots(3, 1,figsize=(12,9))

for i,file in enumerate(files):
	print(file)
	reactionFile=runLocation+file+'reactionHistory.csv'
	# area (mm^2)
	# time, F00, F11, F22, length_x, length_y, length_z, Rx-, Rx+, Ry-, Ry+, Rz-, Rz+, L00, L11, L22
	data = np.genfromtxt(reactionFile, delimiter=',')
	time = data[:,0]
	F00 = data[:,1]
	F11 = data[:,2]
	F22 = data[:,3]
	lenx = data[:,4]
	leny = data[:,5]
	lenz = data[:,6]
	Rxm = data[:,7]
	Rxp = data[:,8]
	Rym = data[:,9]
	Ryp = data[:,10]
	Rzm = data[:,11]
	Rzp = data[:,12]
	Lxx = data[:,13]
	Lyy = data[:,14]
	Lzz = data[:,15]

	sampleHeight=heights[i]
	sampleArea=areas[i]
	
	# this hides all non-monotic time entries, so restart data files are cleaned:
	maxt = 0.0
	mask = np.ones(len(time), dtype=bool)
	for ii,t in enumerate(time):
		if (t<=maxt):
			mask[ii] = False
		else:
			maxt = t
	time = time[mask,...]
	F00 = F00[mask,...]
	F11 = F11[mask,...]
	F22 = F22[mask,...]
	lenx = lenx[mask, ...]
	leny = leny[mask, ...]
	lenz = lenz[mask, ...]
	Rxm = Rxm[mask,...]
	Rxp = Rxp[mask,...]
	Rym = Rym[mask,...]
	Ryp = Ryp[mask,...]
	Rzm = Rzm[mask,...]
	Rzp = Rzp[mask,...]
	Lzz = Lzz[mask,...]

	# strain
	exx=np.log(F00)
	eyy=np.log(F11)
	ezz=np.log(F22)

	sampleHeight = F22[1]
	fzz = 1000 * 0.5*(-Rzp+Rzm)
	fzzp = -1000.0*Rzp
	fzzm = 1000 * Rzm
	disp = (1-F22)*sampleHeight

	# vm=np.sqrt(0.5*( np.square(sxx-syy)+np.square(syy-szz)+np.square(szz-sxx) ) )
	grainStrain = disp/sampleHeight

    # stress is scaled from simulation units (GPa) to plotting units (MPa)
	ax1.plot(time,disp,linestyle='-',color=colors[i],linewidth=1,label=str(labels[i]))

	ax2.plot(grainStrain,szz,linestyle='-',color=colors[i],linewidth=1,label=str(labels[i]))  #kN by default
	
	ax3.plot(disp,fzz,linestyle='-.',color=colors[i],linewidth=1,label=str(labels[i]))
	
# stretch vs time
ax1.set_xlabel(r'time ($\mu$s)', fontsize=16)
ax1.set_ylabel('  ', fontsize=16, labelpad=10)
ax1.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax1.grid()

# Eng Stress vs Eng Strain
ax2.set_xlabel(r'Eng. Strain (disp/grainHeight: $\epsilon_{zz}$)', fontsize=16)
ax2.set_ylabel('  ', fontsize=16, labelpad=10)
# ax2.set_xlim(0,0.0006)
# ax2.set_ylim(0.0,75)
ax2.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax2.grid()

# Force vs Compressive Displacement
ax3.set_xlabel('Compressive displacement (mm)', fontsize=16)
# ax3.set_xlabel(r'Time ($\mu s$)', fontsize=16)
ax3.set_ylabel('  ', fontsize=16, labelpad=10)
ax3.grid()
ax3.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
# ax3.set_xlim(0.0, 0.05)
# ax3.set_ylim(-15, 15)

# Write aligned Y axis labels
fig.text(0.01, 0.0, r'         Compressize force (N)          Eng. Stress ($\sigma_{zz}$, $GPa$)             Disp. ($\vert \, u(t) \, \vert$,  $mm$)', va='bottom', rotation='vertical', fontsize=15)
# fig.text(0.01, 0.0, r'               $fzz^+-fzz^-$          Eng. Stress ($\sigma_{zz}$, $GPa$)             Disp. ($\vert \, u(t) \, \vert$,  $mm$)', va='bottom', rotation='vertical', fontsize=15)
fig.tight_layout()
plt.show()
# fig.savefig("forceVSdispDynamicCylinderVerification.png", bbox_inches="tight")