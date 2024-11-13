# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 09:00:00 2017
@author: homel1
Geometry object functions for the particle file writer.
"""
import numpy as np                   # math stuff
import matplotlib.pyplot as plt

import matplotlib.cm as cm

import matplotlib.pyplot as plt
from cycler import cycler

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


# parent directory of file pathL:
runLocation='/p/lustre1/homel1/geosxRuns/'

files=[
'elasticBlockUni',
]  

labels=[
'elasticBlockUni',
]

styles=[
'-',
'--',
'--',
'-',
'--',
'-',
'--',
'-'
]


# colors=[
# 'blue',
# 'green',
# 'green',
# 'red',
# 'red'
# ]

# input material properties
bulk = 70
shear = 24
C11 = bulk + (4.0/3.0)*shear
EE = 9.0*bulk*shear/(3.0*bulk + shear)
C12 = C11 - 2.0*shear

# initial volume fraction (assumes one material with density rho02)
VF0=np.zeros(len(files))
rho0s=2.7

#initialize plots
evenly_spaced_interval = np.linspace(0, 1, len(files))
colors = [cm.rainbow(x) for x in evenly_spaced_interval]
fig, (ax1, ax2, ax3) = plt.subplots(3, 1,figsize=(15,12))
fig.suptitle('F table benchmark with uniaxial strain & stress')
#fig2, (ax4,ax5,ax6) = plt.subplots(3, 1,figsize=(16,12))

# plot box sum data and compute initial volume fraction
# for i,file in enumerate(files):
# 	print(file)
# 	reactionFile=runLocation+file+'/boxSumHistory.csv'

# 	# time,sxx,syy,szz,sxy,syz,sxz,rho,energy
# 	data = np.genfromtxt(reactionFile, delimiter=',')
# 	time = data[:,0]
# 	sxx = data[:,1]
# 	syy = data[:,2]
# 	szz = data[:,3]
# 	sxy = data[:,4]
# 	syz = data[:,5]
# 	sxz = data[:,6]
# 	rho = data[:,7]
# 	e = data[:,8]

# 	# this hids all non-monotic time entries, so restart data files are cleaned:
# 	maxt = 0.0
# 	mask = np.ones(len(time), dtype=bool)
# 	for ii,t in enumerate(time):
# 		if (t<=maxt):
# 			mask[ii] = False
# 		else:
# 			maxt = t
# 	time = time[mask,...]
# 	sxx = sxx[mask,...]
# 	syy = syy[mask,...]
# 	szz = szz[mask,...]
# 	sxy = sxy[mask,...]
# 	syz = syz[mask,...]
# 	sxz = sxz[mask,...]
# 	rho = rho[mask,...]
# 	e = e[mask,...]


# 	# p=(-1.0/3.0)*(sxx+syy+szz)
# 	# vm=np.sqrt(0.5*( np.square(sxx-syy)+np.square(syy-szz)+np.square(szz-sxx) ) )
# 	ax3.plot(time,1000*sxx,linestyle=':',color=cm.gist_rainbow(0),linewidth=1,label='box: sxx')
# 	ax3.plot(time,1000*syy,linestyle=':',color=cm.gist_rainbow(0.33),linewidth=1,label='box: syy')
# 	ax3.plot(time,1000*szz,linestyle=':',color=cm.gist_rainbow(0.67),linewidth=1,label='box: szz')


# print("VF0 = ",VF0) 


for i,file in enumerate(files):
	rho0=rho0s*VF0[i]
	print(file)
	reactionFile=runLocation+file+'/reactionHistory.csv'
	# time,F00,F11,F22,lx,ly,lz,Rx-,Rx+,Ry-,Ry+,Rz-,Rz+
	data = np.genfromtxt(reactionFile, delimiter=',')
	time = data[:,0]
	F00 = data[:,1]
	F11 = data[:,2]
	F22 = data[:,3]
	lx  = data[:,4]
	ly  = data[:,5]
	lz  = data[:,6]
	Rxm = data[:,7]
	Rxp = data[:,8]
	Rym = data[:,9]
	Ryp = data[:,10]
	Rzm = data[:,11]
	Rzp = data[:,12]

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
	Rxm = Rxm[mask,...]
	Rxp = Rxp[mask,...]
	Rym = Rym[mask,...]
	Ryp = Ryp[mask,...]
	Rzm = Rzm[mask,...]
	Rzp = Rzp[mask,...]

	# strain
	exx=np.log(F00)
	eyy=np.log(F11)
	ezz=np.log(F22)
	
	#stress
	Ax=ly*lz
	Ay=lx*lz
	Az=lx*ly

	sxx=0.5*(Rxp-Rxm)/Ax
	syy=0.5*(Ryp-Rym)/Ay
	szz=0.5*(Rzp-Rzm)/Az

	p=(-1.0/3.0)*(sxx+syy+szz)
	rho=rho0/F00

	vm=np.sqrt(0.5*( np.square(sxx-syy)+np.square(syy-szz)+np.square(szz-sxx) ) )


    # stress is scaled from simulation units (GPa) to plotting units (MPa)
	ax1.plot(time,F00,linestyle='-',color=cm.gist_rainbow(0),linewidth=1,label='F00')
	ax1.plot(time,F11,linestyle='-',color=cm.gist_rainbow(0.33),linewidth=1,label='F11')
	ax1.plot(time,F22,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=1,label='F22')

	ax2.plot(exx,-1000*Rxm/Ax,linestyle='-',color=cm.gist_rainbow(0),linewidth=2,label='x_minus')
	ax2.plot(exx,1000.0*Rxp/Ax,linestyle='--',color=lighten_color(cm.gist_rainbow(0),1.5),linewidth=2,label='x_plus')
	ax2.plot(exx,-1000*Rym/Ay,linestyle='-',color=cm.gist_rainbow(0.33),linewidth=2,label='y_minus')
	ax2.plot(exx,1000.0*Ryp/Ay,linestyle='--',color=lighten_color(cm.gist_rainbow(0.33),1.5),linewidth=2,label='y_plus')
	ax2.plot(exx,-1000*Rzm/Az,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=1,label='z_minus')
	ax2.plot(exx,1000.0*Rzp/Az,linestyle='--',color=lighten_color(cm.gist_rainbow(0.67),1.5),linewidth=2,label='z_plus')

	ax2.plot(exx,1e3*C11*exx,linestyle='-',color=cm.gnuplot2(0),linewidth=1,label='Constrained (C11) Modulus')
	ax2.plot(exx,1e3*EE*exx,linestyle='--',color=cm.gnuplot2(0),linewidth=1,label='Unconstrained (Young\'s) Modulus')
	ax2.plot(exx,1e3*C12*exx,linestyle='-.',color=cm.gnuplot2(0),linewidth=1,label='C12')

	ax3.plot(time,-1000*Rxm/Ax,linestyle='-',color=cm.gist_rainbow(0),linewidth=2,label='x_minus')
	ax3.plot(time,1000.0*Rxp/Ax,linestyle='--',color=lighten_color(cm.gist_rainbow(0),1.5),linewidth=2,label='x_plus')
	ax3.plot(time,-1000*Rym/Ay,linestyle='-',color=cm.gist_rainbow(0.33),linewidth=2,label='y_minus')
	ax3.plot(time,1000.0*Ryp/Ay,linestyle='--',color=lighten_color(cm.gist_rainbow(0.33),1.5),linewidth=2,label='y_plus')
	ax3.plot(time,-1000*Rzm/Az,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=2,label='z_minus')
	ax3.plot(time,1000.0*Rzp/Az,linestyle='--',color=lighten_color(cm.gist_rainbow(0.67),1.5),linewidth=2,label='z_plus')



# von Mises vs density
ax1.set_xlabel('time')
#ax1.set_xlim(0.0,-1.05*min(exx[1:]))
ax1.set_ylabel('F')
#ax1.set_yscale("log")
#ax1.set_ylim(0.001,5.0)
#ax1.set_ylim(-2.5,2.5)
ax1.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='x-small')
ax1.grid()

# plot solid eos
# js = np.arange(0.8, 1.0001, 0.0001)
# ps = -bulk*np.log(js)
# ax2.plot(1/js, ps, 'k-',linewidth=2,label='Solid EOS')

ax2.set_xlabel('Applied Strain (mm/mm)', fontsize=16)
#ax2.set_xlim(0,1.2)
#ax2.set_ylim(0,600)
#ax2.set_yscale("log")
#ax2.set_ylim(0.001,1000)

ax2.set_ylabel(r'$\sigma$ (MPa)', fontsize=16)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax2.tick_params(axis='both', which='minor', labelsize=16)

ax2.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax2.grid()


ax3.set_xlabel('time', fontsize=16)
ax3.set_ylabel(r'$\sigma$ (MPa)', fontsize=16)
ax3.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')

ax3.set_xlim(0.00,25)
ax3.set_ylim(0.00,200)
#ax3.set_yscale("log")
#ax3.set_ylim(0.001,1000)


fig.tight_layout()
#plt.show()
fig.savefig("elasticBlockUni.png", bbox_inches="tight")