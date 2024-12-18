import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
from xml.etree import ElementTree
import dpWellbore_analyticalResults as dpAnal
import elastoPlasticWellboreAnalyticalSolutions as epwAnal

def extractDataFromXMLList(paramList):
	# Extract data from a list in XML such as "{ 1, 2, 3}"
	return paramList.replace('{', '').replace('}', '').strip().split(',')

def getDataFromXML(xmlFilePathPrefix):
	# Get wellbore inner radius
	xmlFilePath = xmlFilePathPrefix + "_benchmark.xml"
	tree = ElementTree.parse(xmlFilePath)

	meshParam = tree.find('Mesh/InternalWellbore')
	radii = extractDataFromXMLList( meshParam.get("radius") )
	a0 = float(radii[0])

	# Get the temperature change on the inner surface of the wellbore
	xmlFilePath = xmlFilePathPrefix + "_base.xml"
	tree = ElementTree.parse(xmlFilePath)

	fsParams = tree.findall('FieldSpecifications/FieldSpecification')
	for fsParam in fsParams:
		if fsParam.get('name') == 'stressXX':
			sh = float( fsParam.get('scale') )
		if fsParam.get('name') == 'stressZZ':
			sv = float( fsParam.get('scale') )

	pw = float( extractDataFromXMLList(tree.find('Functions/TableFunction').get('values') )[1])

	defaultBulkModulus = float( tree.find('Constitutive/DruckerPrager').get('defaultBulkModulus') )
	defaultShearModulus = float( tree.find('Constitutive/DruckerPrager').get('defaultShearModulus') )
	defaultCohesion = float( tree.find('Constitutive/DruckerPrager').get('defaultCohesion') )
	defaultFrictionAngle = float( tree.find('Constitutive/DruckerPrager').get('defaultFrictionAngle') )
	defaultDilationAngle = float( tree.find('Constitutive/DruckerPrager').get('defaultDilationAngle') )
	defaultHardeningRate = float( tree.find('Constitutive/DruckerPrager').get('defaultHardeningRate') )

	return [a0, sh, sv, pw, defaultBulkModulus, defaultShearModulus, defaultCohesion, defaultFrictionAngle, defaultDilationAngle, defaultHardeningRate]

def stressRotation(stress, phi_x):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x), 0.], [-np.sin(phi_x), np.cos(phi_x), 0.], [0., 0., 1.]])
    return np.dot(np.dot(np.transpose(rotx), stress), rotx)

def main():
	xmlFilePathPrefix = "DruckerPragerWellbore"
	xmlData = getDataFromXML(xmlFilePathPrefix)

	# Initial wellbore radius	
	a0 = xmlData[0] 

	# Initial stresses
	sh = -xmlData[1] #negative sign because of positive sign convention for compression stress in the analytical solution
	sv = -xmlData[2] #negative sign because of positive sign convention for compression stress in the analytical solution

	# Boundary condition on the wellbore surface
	pw = -xmlData[3] #negative sign because of positive sign convention for compression stress in the analytical solution 
	# Equivalent displacement condition for the analytical solution
	a0_a_ratio = 1.1 # This ratio corresponds to pw, extracted from the relationship between pw and a0/a (see bottom-right figure) 

	# Elastic moduli
	K = xmlData[4]
	G = xmlData[5]
	nu = (3.0*K-2.0*G) / (6.0*K+2.0*G)    
	
	# Elasto-plastic parameters
	defaultCohesion = xmlData[6]
	defaultFrictionAngle = xmlData[7] # deg
	defaultDilationAngle = xmlData[8]
	defaultHardeningRate = xmlData[9]
	
	# Get stress, time and element center from GEOSX results
	stress_field_name = 'rock_stress'
	hf_stress = h5py.File('stressHistory_rock.hdf5', 'r')
	time = np.asarray( hf_stress.get(stress_field_name + ' Time') )
	center = np.asarray( hf_stress.get(stress_field_name + ' elementCenter') )
	stress = np.asarray( hf_stress.get(stress_field_name) )
	
	# Get the deformed wellbore radius
	hf_disp = h5py.File("displacementHistory.hdf5", 'r')
	displacement = np.asarray( hf_disp.get('totalDisplacement') )
	node_position = np.asarray( hf_disp.get('totalDisplacement ReferencePosition') )
	da = displacement[:, 0, 0]
	a = a0 + da
	
	# Compute total stress, the stress of each element is the average value of its eight Gauss points
	nTimeSteps = 39
	stress_xx = ( sum(stress[nTimeSteps,:,6*i+0] for i in range(8)) )/8.0
	stress_yy = ( sum(stress[nTimeSteps,:,6*i+1] for i in range(8)) )/8.0
	stress_zz = ( sum(stress[nTimeSteps,:,6*i+2] for i in range(8)) )/8.0
	stress_yz = ( sum(stress[nTimeSteps,:,6*i+3] for i in range(8)) )/8.0
	stress_xz = ( sum(stress[nTimeSteps,:,6*i+4] for i in range(8)) )/8.0
	stress_xy = ( sum(stress[nTimeSteps,:,6*i+5] for i in range(8)) )/8.0
	
	# Coordinate of elemnt center
	nElements = center.shape[1]
	xCoord = center[0, :, 0]
	yCoord = center[0, :, 1]

	rCoord = np.sqrt( xCoord*xCoord + yCoord*yCoord )

	# Compute stress components in cylindrical coordinate system
	stress_rr = np.zeros(stress_xx.shape)
	stress_tt = np.zeros(stress_xx.shape)

	for idx_elem in range(stress.shape[1]):
		stressMatrix_cartesian = np.array([[stress_xx[idx_elem],stress_xy[idx_elem],stress_xz[idx_elem]],\
                         				   [stress_xy[idx_elem],stress_yy[idx_elem],stress_yz[idx_elem]],\
                         				   [stress_xz[idx_elem],stress_yz[idx_elem],stress_zz[idx_elem]]])

		if(yCoord[idx_elem] != 0):
			phi_x = np.arctan( xCoord[idx_elem]/yCoord[idx_elem] )
		else:
			phi_x = 0

		stressMatrix_cylindirical = stressRotation(stressMatrix_cartesian, phi_x)
		stress_rr[idx_elem] = stressMatrix_cylindirical[1][1]
		stress_tt[idx_elem] = stressMatrix_cylindirical[0][0]

	# Stress invariants at wellbore surface
	stress_xx_wellboresurface = ( sum(stress[:,0,6*i+0] for i in range(8)) )/8.0
	stress_yy_wellboresurface = ( sum(stress[:,0,6*i+1] for i in range(8)) )/8.0
	stress_zz_wellboresurface = ( sum(stress[:,0,6*i+2] for i in range(8)) )/8.0
	stress_yz_wellboresurface = ( sum(stress[:,0,6*i+3] for i in range(8)) )/8.0
	stress_xz_wellboresurface = ( sum(stress[:,0,6*i+4] for i in range(8)) )/8.0
	stress_xy_wellboresurface = ( sum(stress[:,0,6*i+5] for i in range(8)) )/8.0
	p = (stress_xx_wellboresurface + stress_yy_wellboresurface + stress_zz_wellboresurface)/3.0
	q = np.sqrt(1.5) * np.sqrt( (stress_xx_wellboresurface-p)**2.0 + (stress_yy_wellboresurface-p)**2.0 + (stress_zz_wellboresurface-p)**2.0 + \
        stress_xy_wellboresurface**2.0 + stress_xz_wellboresurface**2.0 + stress_yz_wellboresurface**2.0 )
	
	# Analytical results 
	r_ep_analytic,srr_ep_analytic,stt_ep_analytic,szz_ep_analytic,r_elas_analytic,srr_elas_analytic,stt_elas_analytic,szz_elas_analytic,\
    pw_analytic,p_wellsurface_analytic,q_wellsurface_analytic = dpAnal.DruckerPragerModel(a0_a_ratio, sh, sv, nu, a0, G, defaultCohesion, defaultFrictionAngle, defaultDilationAngle, defaultHardeningRate)
	
	# Compute pressure at wellbore surface
	list_a0_a_ratio = np.arange(1.0,a0_a_ratio,0.001)
	list_pw_analytic = []
	list_p_wellsurface_analytic = []
	list_q_wellsurface_analytic = []
	
	for a0_a_ratio_val in list_a0_a_ratio:
		tmp = dpAnal.DruckerPragerModel(a0_a_ratio_val, sh, sv, nu, a0, G, defaultCohesion, defaultFrictionAngle, defaultDilationAngle, defaultHardeningRate)
		pw_analytic = tmp[8]
		p_wellsurface_analytic = tmp[9]
		q_wellsurface_analytic = tmp[10]
		list_pw_analytic.append(pw_analytic)
		list_p_wellsurface_analytic.append(p_wellsurface_analytic)
		list_q_wellsurface_analytic.append(q_wellsurface_analytic)
	
	# Plots	
	plt.figure(figsize=(15,10))

	# Plot GEOS results versus analytical results for elasto-plastic material
	plt.subplot(2,2,1)
	plt.plot( rCoord/rCoord[0],-stress_rr,'r+',label=r'$\sigma_{rr}$: GEOS')
	plt.plot(r_ep_analytic,srr_ep_analytic, 'r', label=r'$\sigma_{rr}$: Analytical Plasticity')
	plt.plot(r_elas_analytic,srr_elas_analytic, 'r')

	plt.plot( rCoord/rCoord[0],-stress_tt,'b+',label=r'$\sigma_{\theta\theta}$: GEOS')
	plt.plot(r_ep_analytic,stt_ep_analytic, 'b', label=r'$\sigma_{\theta\theta}$: Analytical Plasticity')
	plt.plot(r_elas_analytic,stt_elas_analytic, 'b')

	plt.plot( rCoord/rCoord[0],-stress_zz,'g+',label=r'$\sigma_{zz}$: GEOS')
	plt.plot(r_ep_analytic,szz_ep_analytic, 'g', label=r'$\sigma_{zz}$: Analytical Plasticity')
	plt.plot(r_elas_analytic,szz_elas_analytic, 'g')

	plt.plot([r_ep_analytic[0],r_ep_analytic[0]],[0.0,20e6],'k--', label='Elastic-Plastic boundary')

	plt.ylabel('Stresses [Pa]')
	plt.xlabel(r'Normalized radial coordinate, $r/a$')
	plt.xlim(1.0,10.0)
	plt.ylim(0.0,20e6)
	plt.yticks(np.linspace(0.0,20e6,11))
	plt.xscale('log', subs=range(2,10))
	plt.legend()

	# Plot GEOS results for elasto-plastic material versus analytical results for elastic material
	r_elas_assumed,srr_elas_assumed,stt_elas_assumed,szz_elas_assumed = epwAnal.solution_elastic(sh,sv,1.0,pw)
	
	plt.subplot(2,2,2)
	plt.plot( rCoord/rCoord[0],-stress_rr,'r+',label=r'$\sigma_{rr}$: GEOS')
	plt.plot(r_elas_assumed,srr_elas_assumed, 'r', label=r'$\sigma_{rr}$: Analytical Elasticity')

	plt.plot( rCoord/rCoord[0],-stress_tt,'b+',label=r'$\sigma_{\theta\theta}$: GEOS')
	plt.plot(r_elas_assumed,stt_elas_assumed, 'b', label=r'$\sigma_{\theta\theta}$: Analytical Elasticity')

	plt.plot( rCoord/rCoord[0],-stress_zz,'g+',label=r'$\sigma_{zz}$: GEOS')
	plt.plot(r_elas_assumed,szz_elas_assumed, 'g', label=r'$\sigma_{zz}$: Analytical Elasticity')
	
	plt.ylabel('Stresses [Pa]')
	plt.xlabel(r'Normalized radial coordinate, $r/a$')
	plt.xlim(1.0,10.0)
	plt.ylim(0.0,20e6)
	plt.yticks(np.linspace(0.0,20e6,11))
	plt.xscale('log', subs=range(2,10))
	plt.legend()
	
	# Plot the stress path on the p-q plan when the well surface pressure decrease from the in-situ condition
	defaultFrictionAngle *= np.pi/180.0 # converted to rad
	param_b = dpAnal.compute_param_b(defaultFrictionAngle)
	param_a = defaultCohesion * param_b / np.tan(defaultFrictionAngle)
	p_yieldSurface = np.arange(0,15e6,1)
	q_yieldSurface_i = p_yieldSurface*param_b + param_a

	plt.subplot(2,2,3)
	plt.plot(-p,q,'ko',label='GEOS')
	plt.plot(list_p_wellsurface_analytic,list_q_wellsurface_analytic,'b',label='Analytical')
	plt.plot(p_yieldSurface,q_yieldSurface_i,'g--',label='Initial Yield surface')

	plt.xlabel('p (Pa)')
	plt.ylabel('q (Pa)')
	#plt.xlim(0,15e6)
	#plt.ylim(0,10e6)
	plt.legend()

	# Plot the wellbore deformation a0/a versus wellbore surface pressure pw
	pw = stress_xx_wellboresurface

	plt.subplot(2,2,4)
	plt.plot(a0/a, -pw,'ko',label='GEOS')
	plt.plot(list_a0_a_ratio,list_pw_analytic,'b',label='Analytical')
	plt.xlabel('a0/a')
	plt.ylabel('pw (Pa)')
	plt.xlim(1.0,a0_a_ratio)
	plt.ylim(0,12e6)
	plt.legend()
	
	plt.savefig('dpWellboreVerification.png')

if __name__ == "__main__":
	main()
	
