import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.linalg
from scipy import special
from xml.etree import ElementTree

# Analytical results
def steadyState(Tin, Tout, Rin, Rout, radialCoordinate):
	return Tin + (Tout - Tin) * (np.log(radialCoordinate) - np.log(Rin)) / (np.log(Rout) - np.log(Rin))

def diffusionFunction(radialCoordinate, Rin, diffusionCoefficient, diffusionTime):
	return special.erfc(  (radialCoordinate - Rin) / 2.0 / np.sqrt( diffusionCoefficient * diffusionTime ) )

def computeTransientTemperature(Tin, Rin, radialCoordinate, thermalDiffusionCoefficient, diffusionTime):
	# Ref. Wang and Papamichos (1994), https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94WR01774
	return Tin * np.sqrt(Rin/radialCoordinate) * diffusionFunction(radialCoordinate, Rin, thermalDiffusionCoefficient, diffusionTime)
	 
def computeThermalDiffusionCoefficient(thermalConductivity, volumetricHeatCapacity):
	return thermalConductivity / volumetricHeatCapacity
	
def extractDataFromXMLList(paramList):
	# Extract data from a list in XML such as "{ 1, 2, 3}"
	return paramList.replace('{', '').replace('}', '').strip().split(',')

def getWellboreGeometryFromXML(xmlFilePath):
	tree = ElementTree.parse(xmlFilePath)

	meshParam = tree.find('Mesh/InternalWellbore')
	radii = extractDataFromXMLList( meshParam.get("radius") )

	Rin = float(radii[0])
	Rout = float(radii[-1])

	return [Rin, Rout]

def getLoadingFromXML(xmlFilePath):
	tree = ElementTree.parse(xmlFilePath)
	fsParams = tree.findall('FieldSpecifications/FieldSpecification')

	for fsParam in fsParams:
		if ( (fsParam.get('fieldName') == "pressure") & (fsParam.get('initialCondition') != "1") ): 
			if fsParam.get('setNames') == "{ rneg }":
				Pin = float(fsParam.get('scale'))
			if fsParam.get('setNames') == "{ rpos }":
				Pout = float(fsParam.get('scale'))

	for fsParam in fsParams:
		if ( (fsParam.get('fieldName') == "temperature") & (fsParam.get('initialCondition') != "1") ): 
			if fsParam.get('setNames') == "{ rneg }":
				Tin = float(fsParam.get('scale'))
			if fsParam.get('setNames') == "{ rpos }":
				Tout = float(fsParam.get('scale'))
            
	thermalConductivity = float( extractDataFromXMLList( tree.find('Constitutive/SinglePhaseConstantThermalConductivity').get('thermalConductivityComponents') )[0] )

	tree_SolidInternalEnergies = tree.findall('Constitutive/SolidInternalEnergy')

	for tree_SolidInternalEnergy in tree_SolidInternalEnergies:
		if tree_SolidInternalEnergy.get('name') == "rockInternalEnergy_linear":
			volumetricHeatCapacity = float( tree_SolidInternalEnergy.get('volumetricHeatCapacity') )
			
	
	permeability = float( extractDataFromXMLList( tree.find('Constitutive/ConstantPermeability').get('permeabilityComponents') )[0] )

	porosity = float( tree.find('Constitutive/PressurePorosity').get('defaultReferencePorosity') )
	
	fluidViscosity = float( tree.find('Constitutive/ThermalCompressibleSinglePhaseFluid').get('defaultViscosity') )

	fluidCompressibility = float( tree.find('Constitutive/ThermalCompressibleSinglePhaseFluid').get('compressibility') )

	fluidThermalExpansionCoefficient = float( tree.find('Constitutive/ThermalCompressibleSinglePhaseFluid').get('thermalExpansionCoeff') )

	return [Pin, Pout, Tin, Tout, thermalConductivity, volumetricHeatCapacity, permeability, porosity, fluidViscosity, fluidCompressibility, fluidThermalExpansionCoefficient]


def main():

	xmlFilePath = "../../../../../../../inputFiles/singlePhaseFlow/thermalCompressible_2d"
	
	Rin, Rout = getWellboreGeometryFromXML(xmlFilePath+"_benchmark.xml")

	Pin, Pout, Tin, Tout, thermalConductivity, volumetricHeatCapacity, permeability, porosity, fluidViscosity, fluidCompressibility, fluidThermalExpansionCoefficient = getLoadingFromXML(xmlFilePath+"_base.xml")

	plt.figure(figsize=(10,7))
	font = {'size'   : 16}
	plt.rc('font', **font)

	for chart_idx, idx in enumerate([1, 2, 5, 10]):
		# Numerical results
		data = pd.read_csv(f'data_{idx}.csv')
		radialCoordinate = data['elementCenter:0']
		temperature = data['temperature']
		pressure = data['pressure']
		diffusionTime = data['Time'][0]

		# Analytical results
		thermalDiffusionCoefficient = computeThermalDiffusionCoefficient(thermalConductivity, volumetricHeatCapacity)

		T_transient = computeTransientTemperature(Tin, Rin, radialCoordinate, thermalDiffusionCoefficient, diffusionTime)
		
		# Analytical results of the steady state regime for comparison
		T_steadyState = steadyState(Tin, Tout, Rin, Rout, radialCoordinate)

		# Visualization
		# Temperature
		plt.subplot(2,2,chart_idx+1)
		plt.plot( radialCoordinate, temperature, 'k+' , label='GEOSX' )
		plt.plot( radialCoordinate, T_transient, 'r-' , label='Analytic' )
		plt.plot( radialCoordinate, T_steadyState, 'b-' , label='Steady State' )

		if chart_idx==1:
			plt.legend()

		if chart_idx in [2,3]:
			plt.xlabel('Radial distance from well center')

		if chart_idx in [0,2]:
			plt.ylabel('Temperature (°C)')

		plt.ylim(-10,100)
		plt.xlim(0,1.0)
		plt.title('t = '+str(diffusionTime)+'(s)')
		plt.tight_layout()

	plt.show()

if __name__ == "__main__":
	main()
