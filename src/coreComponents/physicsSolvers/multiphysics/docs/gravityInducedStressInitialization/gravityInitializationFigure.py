import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
import math
from math import sin,cos,tan,exp,atan,asin
import csv


def getHydromechanicalParametersFromXML(xmlFilePath):
	tree = ElementTree.parse(xmlFilePath)

	param1 = tree.find('Constitutive/ElasticIsotropic')
	param2 = tree.find('Constitutive/BiotPorosity')
	param3 = tree.find('Constitutive/CompressibleSinglePhaseFluid')

	hydromechanicalParameters = dict.fromkeys([
            "bulkModulus", "shearModulus", "youngModulus", "poissonRatio", "rockDensity", "poissonRatio", "biotCoefficient", "porosity", "fluidDensity", "traction"])

	hydromechanicalParameters["rockDensity"] = float(param1.get("defaultDensity"))
	hydromechanicalParameters["poissonRatio"] = float(param1.get("defaultPoissonRatio"))
	hydromechanicalParameters["youngModulus"] = float(param1.get("defaultYoungModulus"))

	E = hydromechanicalParameters["youngModulus"] 
	nu = hydromechanicalParameters["poissonRatio"]
	K = E / (3 * (1 - 2 * nu))
	G = E / (2 * (1 + nu))

	hydromechanicalParameters["poissonRatio"] = nu
	hydromechanicalParameters["bulkModulus"] = K
	hydromechanicalParameters["shearModulus"] = G

	Ks = float(param2.get("defaultGrainBulkModulus"))
	hydromechanicalParameters["biotCoefficient"] = 1.0 - K / Ks

	hydromechanicalParameters["porosity"] = float(param2.get("defaultReferencePorosity"))

	hydromechanicalParameters["fluidDensity"] = float(param3.get("defaultDensity"))

	param4 = tree.findall('FieldSpecifications/Traction')
	found_stress = False
	for elem in param4:
		if elem.get("name") == "tractionTop" and elem.get("tractionType") == "normal":
			traction = float(elem.get("scale")) * (-1)
			found_stress = True
		if found_stress: break 

	return hydromechanicalParameters


def main():
	# File path
	xmlFile1Path = "../../../../../../inputFiles/initialization/gravityInducedStress_initialization_base.xml"
	xmlFile2Path = "../../../../../../inputFiles/initialization/gravityInducedStress_initialization_benchmark.xml"

	hydromechanicalParameters = getHydromechanicalParametersFromXML(xmlFile1Path)

	BiotCoefficient = hydromechanicalParameters["biotCoefficient"]
	nu = hydromechanicalParameters["poissonRatio"]
	rhoF = hydromechanicalParameters["fluidDensity"]
	rhoR = hydromechanicalParameters["rockDensity"]
	phi = hydromechanicalParameters["porosity"]
	rhoB = (1-phi)*rhoR + phi*rhoF
	
	traction = hydromechanicalParameters["traction"]
	gravity = 9.81 
	
	# rename this file to the name of your Paraview output file
	file = open("simulation_result_0.csv")
	csvreader = csv.reader(file)
	header = next(csvreader)
	header_index = {column_name: index for index, column_name in enumerate(header)}

	rows = []
	for row in csvreader:
		rows.append(row)
	file.close() 

	zloc_index = header_index["elementCenter:2"]
	pressure_index = header_index["pressure"]
	tsxx_index = header_index["rockSolid_stress:0"] # the solidModelName="rockSolid" has been defined in the gravityInducedStress_initialization_base.xml file, please change if you have a different solidModelName 
	tsyy_index = header_index["rockSolid_stress:1"]
	tszz_index = header_index["rockSolid_stress:2"]
	

	rows = np.array(rows)
	zloc_0 = np.empty(len(rows[:,1]))
	pressure_0 = np.empty(len(rows[:,1]))
	tsxx_0 = np.empty(len(rows[:,1]))
	tsyy_0 = np.empty(len(rows[:,1]))
	tszz_0 = np.empty(len(rows[:,1]))
	for i in range(0,len(rows[:,1])):
		zloc_0[i]=-(float(rows[i,zloc_index]))
		pressure_0[i]=float(rows[i,pressure_index])		
		tsxx_0[i]=-(float(rows[i,tsxx_index])-BiotCoefficient*pressure_0[i])/1.0e6
		tsyy_0[i]=-(float(rows[i,tsyy_index])-BiotCoefficient*pressure_0[i])/1.0e6
		tszz_0[i]=-(float(rows[i,tszz_index])-BiotCoefficient*pressure_0[i])/1.0e6


	z_analytical= np.linspace(0, 1000, 100)
	pp_analytical= rhoF*gravity*z_analytical/1.0e6
	szz_analtyical= rhoB*gravity*z_analytical/1.0e6

	sxx_analtyical=nu/(1-nu)*(szz_analtyical-BiotCoefficient*pp_analytical)+BiotCoefficient*pp_analytical

	fsize = 20
	msize = 12
	lw = 6
	mew = 2
	malpha = 0.6
	lalpha = 0.8
	N1=1

	fig = plt.figure(figsize=(10,8))
	cmap = plt.get_cmap("tab10")

		
	plt.plot(tsxx_0[::N1], zloc_0[::N1], 'o', color=cmap(0), markersize=msize, alpha=malpha, mec=cmap(0), fillstyle='none', mew=mew, label= 'Sxx_Total_GEOS')
	plt.plot(sxx_analtyical, z_analytical, lw=lw, alpha=0.8, color='orange', linestyle= ':', label='Sxx_Total_Analytical')
	plt.plot(tsyy_0[::N1], zloc_0[::N1], 's', color=cmap(1), markersize=msize, alpha=malpha, mec=cmap(1), fillstyle='none', mew=mew, label= 'Syy_Total_GEOS')
	plt.plot(tszz_0[::N1], zloc_0[::N1], 'd', color=cmap(2), markersize=msize, alpha=malpha, mec=cmap(2), fillstyle='none', mew=mew, label= 'Szz_Total_GEOS')
	plt.plot(szz_analtyical, z_analytical, lw=lw, alpha=0.8, color='g', linestyle= ':', label='Szz_Total_Analytical')
	plt.plot(pressure_0[::N1]/1.0e6, zloc_0[::N1], 'x', color=cmap(3), markersize=msize, alpha=malpha, mec=cmap(3), fillstyle='none', mew=mew, label= 'Pore Pressure_GEOS')
	plt.plot(pp_analytical, z_analytical, lw=lw, alpha=0.8, color='r', linestyle= ':', label='Pore Pressure_Analytical')
	plt.xlabel('Total Stresses [MPa]', size=fsize, weight="bold")
	plt.ylabel('Depth [m]', size=fsize, weight="bold")
	plt.legend(loc='upper right',fontsize=fsize*0.5)
	plt.grid(True)
	ax = plt.gca()
	ax.xaxis.set_tick_params(labelsize=fsize)
	ax.yaxis.set_tick_params(labelsize=fsize)
	ax.invert_yaxis()
	
	plt.show() 


if __name__ == "__main__":
		main()




		


