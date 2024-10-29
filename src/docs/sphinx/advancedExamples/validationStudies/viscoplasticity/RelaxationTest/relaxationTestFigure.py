import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
from math import sin,cos,tan,exp,atan
import h5py


def getMechanicalParametersFromXML( xmlFilePath ):
	tree = ElementTree.parse(xmlFilePath)
    
	param = tree.find('Constitutive/ViscoExtendedDruckerPrager')

	mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus", "cohesion", "initialFrictionAngle", "residualFrictionAngle", "dilationRatio", "hardening", "relaxationTime"])
	mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
	mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
	mechanicalParameters["cohesion"] = float(param.get("defaultCohesion"))
	mechanicalParameters["initialFrictionAngle"] = float(param.get("defaultInitialFrictionAngle"))
	mechanicalParameters["residualFrictionAngle"] = float(param.get("defaultResidualFrictionAngle"))
	mechanicalParameters["dilationRatio"] = float(param.get("defaultDilationRatio"))
	mechanicalParameters["hardening"] = float(param.get("defaultHardening"))
	mechanicalParameters["relaxationTime"] = float(param.get("relaxationTime"))

	return mechanicalParameters


def getLoadingFromXML(xmlFilePath):
	tree = ElementTree.parse(xmlFilePath)

	param = tree.findall('FieldSpecifications/FieldSpecification')
	
	for elem in param:
		if elem.get("name") == "stressZZ" and elem.get("initialCondition") == "1":
			initialStress = float(elem.get("scale"))

		elif elem.get("name") == "axialload" and elem.get("fieldName") == "totalDisplacement":
			strainScale = float(elem.get("scale"))

	param1 = tree.find('Functions/TableFunction')
	if param1.get("name") == "timeFunction" and param1.get("inputVarNames") == "{ time }":
		load = param1.get("values")
		load = [float(i)*strainScale for i in load[1:-1].split(",")]		
		loadTime = param1.get("coordinates")
		loadTime = [float(i) for i in loadTime[1:-1].split(",")]

	return initialStress, load, loadTime  


def semiAnalytical( mechanicalParameters, imp_strain, imp_time, initialStress ):
	bulkModulus = mechanicalParameters["bulkModulus"]
	shearModulus = mechanicalParameters["shearModulus"]
	cohesion = mechanicalParameters["cohesion"]
	dilationRatio = mechanicalParameters["dilationRatio"]
	hardeningParameter = mechanicalParameters["hardening"]
	relaxationTime = mechanicalParameters["relaxationTime"]
	# Compute Lame modulus and Young modulus
	lameModulus = bulkModulus - 2.0/3.0*shearModulus
	youngModulus = 1.0/(1.0/9.0/bulkModulus + 1.0/3.0/shearModulus)

	# Initial friction and cohesion parameters
	initialFrictionAngleRad = mechanicalParameters["initialFrictionAngle"]*np.pi/180.0
	cosInitialFrictionAngle = np.cos(initialFrictionAngleRad)
	sinInitialFrictionAngle = np.sin(initialFrictionAngleRad) 
	a_init = 6.0*cohesion*cosInitialFrictionAngle/(3.0-sinInitialFrictionAngle)
	b_init = 6.0*sinInitialFrictionAngle/(3.0-sinInitialFrictionAngle)

	# Residual friction parameters
	residualFrictionAngleRad = mechanicalParameters["residualFrictionAngle"]*np.pi/180.0
	sinResidualFrictionAngle = np.sin(residualFrictionAngleRad) 
	b_resi = 6.0*sinResidualFrictionAngle/(3.0-sinResidualFrictionAngle)

	list_ax_strain_anal = []
	numStepPerLoadingPeriod = 1000

	for i in range(0,len(imp_strain)-1):
		loadingPeriod = np.linspace(imp_strain[i],imp_strain[i+1], numStepPerLoadingPeriod, endpoint=False)
		list_ax_strain_anal = np.concatenate((list_ax_strain_anal, loadingPeriod), axis=0)

	list_time_anal = []
	for i in range(0,len(imp_time)-1):
		timePeriod = np.linspace(imp_time[i],imp_time[i+1], numStepPerLoadingPeriod, endpoint=False)
		list_time_anal = np.concatenate((list_time_anal, timePeriod), axis=0)

	list_ra_stress_anal = initialStress*np.ones(len(list_ax_strain_anal)) #constant radial confining stress

	# Initiate radial strain and axial stress arrays
	list_ra_strain_anal = np.zeros(len(list_ax_strain_anal))
	list_ax_stress_anal = np.zeros(len(list_ax_strain_anal))
	list_ax_stress_anal[0] = initialStress

	# Loop over the loading/unloading steps	
	list_ra_strain_anal[0] = 0
	plasticMultiplier = 0
	b = b_init

	for idx in range(1,len(list_ax_strain_anal)):
		delta_ax_strain_anal = list_ax_strain_anal[idx] - list_ax_strain_anal[idx-1]
		delta_time_anal = list_time_anal[idx]-list_time_anal[idx-1]
		delta_ra_stress_anal = 0

		# Elastic trial
		delta_ra_strain_anal = (delta_ra_stress_anal-lameModulus*delta_ax_strain_anal)/(2.0*lameModulus+2.0*shearModulus)
		delta_ax_stress_anal = (lameModulus+2.0*shearModulus)*delta_ax_strain_anal + lameModulus/(lameModulus+shearModulus)*(delta_ra_stress_anal-lameModulus*delta_ax_strain_anal)

		ax_stress_anal = list_ax_stress_anal[idx-1] + delta_ax_stress_anal
		ra_strain_anal = list_ra_strain_anal[idx-1] + delta_ra_strain_anal

		# Compute mean and shear stresses	
		ra_stress_anal = list_ra_stress_anal[idx]
		p_anal = (ax_stress_anal + 2.0 * ra_stress_anal) / 3.0
		q_anal = -(ax_stress_anal - ra_stress_anal)

		# Plastic correction
		if(q_anal>=0): #loading
			
			F_anal = q_anal + b*p_anal - b*a_init/b_init

			if(F_anal>=0):
				
				b = b_init + plasticMultiplier/(hardeningParameter+plasticMultiplier) * (b_resi - b_init)
				b_dilation = b*dilationRatio

				dF_db = p_anal - a_init/b_init
				db_dlambda = hardeningParameter * (b_resi - b_init) / (hardeningParameter+plasticMultiplier) / (hardeningParameter+plasticMultiplier)
				hardeningRate = -dF_db*db_dlambda

				# Variation of Perzyna plastic multiplier, see Runesson et al. 1999, see Eq. 56, 4, 80, 62, 63
				parameter_Aep = 3.0*shearModulus + bulkModulus*b*b_dilation + hardeningRate
				delta_lambda = delta_time_anal / relaxationTime * (F_anal/parameter_Aep)
				
				# Compute stress and strain variations
				delta_ax_stress_anal = ( delta_ax_strain_anal-delta_lambda*(b_dilation-3.0)/3.0 ) * youngModulus 
				delta_ra_strain_anal = delta_ax_strain_anal	 - 	delta_ax_stress_anal / 2.0 / shearModulus + 3.0/2.0*delta_lambda

				# Compute stress and strain at actual loading step
				ax_stress_anal = list_ax_stress_anal[idx-1] + delta_ax_stress_anal
				ra_strain_anal = list_ra_strain_anal[idx-1] + delta_ra_strain_anal

				# Update plastic multiplier		
				plasticMultiplier += delta_lambda
		
		else: #unloading
			
			F_anal = -q_anal + b*p_anal - b*a_init/b_init # negative sign added to q for q<0 to obtain the absolute value

			if(F_anal>=0):
				b = b_init + plasticMultiplier/(hardeningParameter+plasticMultiplier) * (b_resi - b_init)
				b_dilation = b*dilationRatio

				dF_db = p_anal - a_init/b_init
				db_dlambda = hardeningParameter * (b_resi - b_init) / (hardeningParameter+plasticMultiplier) / (hardeningParameter+plasticMultiplier)
				hardeningRate = -dF_db*db_dlambda
				
				# Variation of Perzyna plastic multiplier, see Runesson et al. 1999, see Eq. 56, 4, 80, 62, 63
				parameter_Aep = 3.0*shearModulus + bulkModulus*b*b_dilation + hardeningRate
				delta_lambda = delta_time_anal / relaxationTime * (F_anal/parameter_Aep)

				# Compute stress and strain variations
				delta_ax_stress_anal = ( delta_ax_strain_anal-delta_lambda*(b_dilation+3.0)/3.0 ) * youngModulus 
				delta_ra_strain_anal = delta_ax_strain_anal	 - 	delta_ax_stress_anal / 2.0 / shearModulus - 3.0/2.0*delta_lambda

				# Compute stress and strain at actual loading step
				ax_stress_anal = list_ax_stress_anal[idx-1] + delta_ax_stress_anal
				ra_strain_anal = list_ra_strain_anal[idx-1] + delta_ra_strain_anal
		
				# Update plastic multiplier		
				plasticMultiplier += delta_lambda	

		list_ax_stress_anal[idx] = ax_stress_anal
		list_ra_strain_anal[idx] = ra_strain_anal
	
        # Preparing data for visualizing semi-analytical results	
	list_p_anal = (list_ax_stress_anal + 2.0 * list_ra_stress_anal) / 3.0
	list_q_anal = -(list_ax_stress_anal - list_ra_stress_anal)

	return list_ax_stress_anal,list_time_anal,list_ax_strain_anal,list_ra_strain_anal,list_p_anal,list_q_anal


def main():
	# File path
	hdf5File1Path = "displacement_history.hdf5"
	hdf5File2Path = "stress_history.hdf5"
	xmlFile1Path = "../../../../../../../inputFiles/solidMechanics/viscoExtendedDruckerPrager_relaxation_base.xml"
	xmlFile2Path = "../../../../../../../inputFiles/solidMechanics/viscoExtendedDruckerPrager_relaxation_benchmark.xml"

	# Read HDF5
	# Global Coordinate of Nodal Point
	hf = h5py.File(hdf5File1Path, 'r')
	xl_node = hf.get('totalDisplacement ReferencePosition topPoint')        
	xcord_node = xl_node[0,:,0]
	ycord_node = xl_node[0,:,1]
	zcord_node = xl_node[0,:,2]
	# Load Displacement Components
	disp = hf.get('totalDisplacement topPoint')    

	# Global Coordinate of Element Center
	hf = h5py.File(hdf5File2Path, 'r')
	xl_elm = hf.get('rock_stress elementCenter')
	xl_elm = np.array(xl_elm)        
	xcord_elm = xl_elm[0, :, 0]
	ycord_elm = xl_elm[0, :, 1]
	zcord_elm = xl_elm[0, :, 2]
	time = hf.get('rock_stress Time')
	# Load Stress Components
	sigma = hf.get('rock_stress')
	sigma = np.array(sigma)
	sigma_Cart = np.zeros([len(sigma[:,0,0]), len(sigma[0,:,0]), 6]) 
	for tl in range(0,len(sigma[:,0,0])):       
		for i in range(0,len(sigma[0,:,0])):
			for j in range(0, 6):
				for k in range(0, 8):
					sigma_Cart[tl,i,j] += sigma[tl, i, j+6*k]/8.

	ax_stress = sigma_Cart[:, -1, 2]
	ra_stress = sigma_Cart[:, -1, 0] 
	p_num = (ax_stress + 2.0 * ra_stress) / 3.0
	q_num = -(ax_stress - ra_stress)
	ax_strain = disp[:, -1, 2] 
	ra_strain = disp[:, -1, 0]


	# Extract info from XML
	mechanicalParameters = getMechanicalParametersFromXML(xmlFile1Path)
	initialStress, imp_strain, imp_time = getLoadingFromXML(xmlFile1Path)
	list_ax_stress_anal,list_time_anal,list_ax_strain_anal,list_ra_strain_anal,list_p_anal,list_q_anal = semiAnalytical( mechanicalParameters, imp_strain, imp_time, initialStress )


	#Visualization parameters
	fsize = 20
	msize = 12
	lw = 6
	malpha = 0.5
	N1 = 3
	fig, ax = plt.subplots(3,1,figsize=(10, 18))
	cmap = plt.get_cmap("tab10")


	# Plot strain versus shear stress
	ax[0].plot(-ax_strain[::N1]*100, q_num[::N1]*1e-6, 'o', color=cmap(2), alpha=malpha, fillstyle='full', markersize=msize, label='GEOS')
	ax[0].plot(-ra_strain[::N1]*100, q_num[::N1]*1e-6, 'o', color=cmap(2), alpha=malpha, fillstyle='full', markersize=msize)
	ax[0].plot(-list_ax_strain_anal*100, list_q_anal*1e-6, lw=lw, alpha=0.8, color=cmap(2), label='Analytical')
	ax[0].plot(-list_ra_strain_anal*100, list_q_anal*1e-6, lw=lw, alpha=0.8, color=cmap(2))
	ax[0].set_xlim([-0.12, 0.12])
	ax[0].set_ylim(0, 12)
	ax[0].set_xlabel(r'Strain (%)', size=fsize, weight="bold")
	ax[0].set_ylabel(r'Deviatoric Stress (MPa)', size=fsize, weight="bold")
	ax[0].legend(loc='lower left',fontsize=fsize)
	ax[0].grid(True, which="both", ls="-")
	ax[0].xaxis.set_tick_params(labelsize=fsize)
	ax[0].yaxis.set_tick_params(labelsize=fsize)


	# Plot stress path and yield surfaces
	phi_i = mechanicalParameters["initialFrictionAngle"]
	phi_r= mechanicalParameters["residualFrictionAngle"]
	c_i = mechanicalParameters["cohesion"]/1.0e6
	f_i = atan(6.0*sin(phi_i/180*np.pi)/(3.0-sin(phi_i/180*np.pi)))*180/np.pi
	f_r = atan(6.0*sin(phi_r/180*np.pi)/(3.0-sin(phi_r/180*np.pi)))*180/np.pi
	d_i = 6.0*c_i*cos(phi_i/180*np.pi)/(3.0-sin(phi_i/180*np.pi))
	po = d_i/tan(f_i/180*np.pi)
	d_r = po*tan(f_r/180*np.pi)

	k_i = tan(f_i/180*np.pi)
	p_Yield = np.linspace(0, 50, 100)
	q_iniYield = k_i*p_Yield+d_i
	k_r = tan(f_r/180*np.pi)
	q_resYield = k_r*p_Yield+d_r

	ax[1].plot(-p_num[::N1]*1e-6, q_num[::N1]*1e-6, 'o', color=cmap(2), alpha=malpha, fillstyle='full', markersize=msize, label='GEOS')
	ax[1].plot(-list_p_anal*1e-6, list_q_anal*1e-6, lw=lw, alpha=0.8, color=cmap(2), label='Analytical')
	ax[1].plot(p_Yield, q_iniYield, lw=lw, alpha=0.8, color='k', linestyle= '--', label='Initial Yield Surface')
	ax[1].plot(p_Yield, q_resYield, lw=lw, alpha=0.8, color='orange', linestyle= '--', label='Residual Yield Surface')
	ax[1].set_xlim([0, 15])
	ax[1].set_ylim([0, 15])
	ax[1].set_xlabel(r'p (MPa)', size=fsize, weight="bold")
	ax[1].set_ylabel(r'q (MPa)', size=fsize, weight="bold")
	ax[1].legend(loc='upper left',fontsize=fsize)
	ax[1].grid(True, which="both", ls="-")
	ax[1].xaxis.set_tick_params(labelsize=fsize)
	ax[1].yaxis.set_tick_params(labelsize=fsize)



	# Plot axial stress and strain with time
	ax[2].plot(time[::N1]/86400, -ax_strain[::N1]*100, 'o', color=cmap(1), alpha=malpha, fillstyle='full', markersize=msize, label='Axial Strain_GEOS')
	ax[2].plot(list_time_anal/86400, -list_ax_strain_anal*100,  lw=lw, alpha=0.8, color=cmap(1), label='Axial Strain_Analytical')
	ax[2].set_xlabel(r'Time (D)', size=fsize, weight="bold")
	ax[2].set_ylabel(r'Axial Strain (%)', size=fsize, weight="bold")
	ax[2].set_ylim([0, 0.2])
	ax[2].legend(loc='upper left',fontsize=fsize)
	ax[2].grid(True, which="both", ls="-")
	ax[2].xaxis.set_tick_params(labelsize=fsize)
	ax[2].yaxis.set_tick_params(labelsize=fsize)
	ax2=ax[2].twinx()
	ax2.plot(time[::N1]/86400, -ax_stress[::N1]/1.0e6, 'o', color=cmap(0), alpha=malpha, fillstyle='full', markersize=msize, label='Axial Stress_GEOS')
	ax2.plot(list_time_anal/86400, -list_ax_stress_anal/1.0e6,  lw=lw, alpha=0.8, color=cmap(0), label='Axial Stress_Analytical')
	ax2.set_ylabel(r'Stress (MPa)', size=fsize, weight="bold")
	ax2.set_ylim([10, 25])
	ax2.legend(loc='lower right',fontsize=fsize)
	ax2.yaxis.set_tick_params(labelsize=fsize)

	

	plt.subplots_adjust(left=0.15, bottom=0.1, right=0.85, top=0.95, wspace=0.4, hspace=0.4)
	plt.show()
	
if __name__ == "__main__":
    main()
