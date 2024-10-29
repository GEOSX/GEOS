import sys
sys.path.append('../')
import numpy as np
import h5py
from xml.etree import ElementTree

import wellboreAnalyticalSolutions

def extractDataFromXMLList(paramList):
	# Extract data from a list in XML such as "{ 1, 2, 3}"
	return paramList.replace('{', '').replace('}', '').strip().split(',')

def getDataFromXML(xmlFilePathPrefix):
	# Get wellbore inner radius
	xmlFilePath = xmlFilePathPrefix + "_benchmark.xml"
	tree = ElementTree.parse(xmlFilePath)

	meshParam = tree.find('Mesh/InternalWellbore')
	radii = extractDataFromXMLList( meshParam.get("radius") )
	ri = float(radii[0])

	# Get the temperature change on the inner surface of the wellbore
	xmlFilePath = xmlFilePathPrefix + "_base.xml"
	tree = ElementTree.parse(xmlFilePath)

	fsParams = tree.findall('FieldSpecifications/FieldSpecification')

	for fsParam in fsParams:
		if ( (fsParam.get('fieldName') == "temperature") & (fsParam.get('initialCondition') != "1") ): 
			if fsParam.get('setNames') == "{ rneg }":
				Ti = float(fsParam.get('scale'))
	
	drainedBulkModulusRock = float( tree.find('Constitutive/ElasticIsotropic').get('defaultBulkModulus') )
	defaultShearModulus = float( tree.find('Constitutive/ElasticIsotropic').get('defaultShearModulus') )
	defaultDrainedLinearTEC = float( tree.find('Constitutive/ElasticIsotropic').get('defaultDrainedLinearTEC') )
	defaultReferencePorosity = float( tree.find('Constitutive/BiotPorosity').get('defaultReferencePorosity') )		
	grainBulkModulus = float( tree.find('Constitutive/BiotPorosity').get('grainBulkModulus') )		
	fluidCompressibility = float( tree.find('Constitutive/ThermalCompressibleSinglePhaseFluid').get('compressibility') )
	fluidViscosity = float( tree.find('Constitutive/ThermalCompressibleSinglePhaseFluid').get('defaultViscosity') )
	fluidThermalExpansionCoefficient = float( tree.find('Constitutive/ThermalCompressibleSinglePhaseFluid').get('thermalExpansionCoeff') )
	thermalConductivity = float( extractDataFromXMLList( tree.find('Constitutive/SinglePhaseConstantThermalConductivity').get('thermalConductivityComponents') )[0] )
	volumetricHeatCapacity = float( tree.find('Constitutive/SolidInternalEnergy').get('volumetricHeatCapacity') )	
	permeability = float( extractDataFromXMLList( tree.find('Constitutive/ConstantPermeability').get('permeabilityComponents') )[0] )

	return [ri, Ti, drainedBulkModulusRock, defaultShearModulus, defaultDrainedLinearTEC, defaultReferencePorosity, grainBulkModulus, fluidCompressibility, fluidViscosity, fluidThermalExpansionCoefficient, permeability, thermalConductivity, volumetricHeatCapacity]

def analyticalResults(t):
	xmlFilePathPrefix = "../../../../../../../inputFiles/wellbore/ThermoPoroElasticWellbore"

	xmlData = getDataFromXML(xmlFilePathPrefix)

	# Geometry and loading
	ri = xmlData[0] # inner radius of the wellbore
	Ti = xmlData[1]	# temperature change on the inner surface of the wellbore

	# Material properties
	# Notations and data by Cheng (2016)
	# Data of Rock salt

	K = xmlData[2] # drained bulk modulus
	G = xmlData[3] # drained shear modulus
	beta_d = 3 * xmlData[4] # drained volumetric thermal expansion coefficient of the porous medium, 3 times of the drained linear thermal expansion coefficient.
	porosity = xmlData[5]
	Ks = xmlData[6]	# bulk modulus of the solid skeleton

	cf = xmlData[7] # fluid compressibility
	muf = xmlData[8] # fluid viscosity
	beta_f = xmlData[9] # volumetric thermal expansion coefficient of fluid   

	permeability = xmlData[10] # intrinsic permeability of rock
	kT = xmlData[11] # thermal conductivity of rock
	volumetricHeatCapacity = xmlData[12] # product of specific heat and density of rock
	
	# Compute parameters for the analytical solutions
	nu = (3.0*K-2.0*G)/(6.0*K+2.0*G)
	E = 2. * G * (1. + nu)
	alpha = 1.0 - K / Ks
	M = 1.0/( porosity*cf + (alpha-porosity)/Ks )
	M11 = K + 4. / 3. * G
	kappa = permeability/muf
	kappaT = kT/volumetricHeatCapacity
	c =  kappa * (M * M11 / (M11 + alpha**2. * M))	
	alpha_d = K*beta_d

	Ku = K + M*alpha*alpha
	S = (3.0*Ku + 4.0*G) /M /(3.0*K+4.0*G)

	beta_s = beta_d # TODO: update for the case porosityTEC != drainedLinearTEC
	beta_v = porosity*(beta_f - beta_s)
	beta_e = beta_d*alpha + beta_v
	
	beta_c = beta_e - 3.0*alpha *K /(3.0*K+4.0*G) *beta_d
	alpha_e = beta_c/S
	
	# Compute analytical results
	r = np.arange(ri, 10. * ri, 0.01 * ri)

	T_thermal, p_thermal, ur_thermal, sig_rr_thermal, sig_tt_thermal, sig_zz_thermal = wellboreAnalyticalSolutions.inTime_thermal(t, r, ri, Ti, G, nu, alpha, alpha_d, c, kappaT, alpha_e)


	return [ r, T_thermal, p_thermal, ur_thermal, sig_rr_thermal, sig_tt_thermal, sig_zz_thermal ]

