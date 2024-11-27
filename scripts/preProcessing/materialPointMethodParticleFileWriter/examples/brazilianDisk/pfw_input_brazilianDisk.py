# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree
import math

pfw = {} 
pfw["runDebug"] = True
stopTime = 2.0

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.648
bulk = 36.3
shear = 26.0
tensileStrength = 0.449
compressiveStrength = 2.27
maximumStrength = 5.0
crackSpeed = 1.8
damagedMaterialFrictionSlope = 0.5773502691896258

# Domain ---------------------------------------------------------------------------------

refine=6
cpp=14
pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=round(pfw["xpar"]*cpp*1.5)   	# grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp 	# grid cells in the y-direction
pfw["nK"]=3  			# grid cells in the z-direction
pfw["ppc"]=2   		# particles per cell in each direction

domainHeight = 1.0
domainWidth = 1.5*domainHeight
domainLength = domainHeight*(pfw["nK"])/(pfw["nJ"])  # m, to get cubic cells

# Define all the geometric objects -------------------------------------------------------
pfw["xmin"] =-0.5*domainWidth	# m
pfw["xmax"] = 0.5*domainWidth	# m
pfw["ymin"] = 0.0	# m
pfw["ymax"] = domainHeight	# m
pfw["zmin"] =-0.5*domainLength 	# m
pfw["zmax"] = 0.5*domainLength 	# m

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime            
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20 # Don't need restarts for now

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25 
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=1
pfw["planeStrain"] = 1

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=1
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["useEvents"]=1
pfw["frictionCoefficient"]=0.25

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Deformation ---------------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 1, 1 ]

pfw["fTableInterpType"]=2
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,        1.00, 1.00, 1.00],
		       [stopTime, 1.00, 0.80, 1.00]]

# Define all the geometric objects -------------------------------------------------------

disk1 = geom.cylinder('disk1',[0,domainHeight/2,pfw["zmin"]],[0,domainHeight/2,pfw["zmax"]],domainHeight/2,[0,0,0],0,0,0)
pfw["objects"]=[disk1]

pfw["materials"] = [ "sand" ]
pfw["materialPropertyString"]="""
<CeramicDamage
	name="sand"
	defaultDensity="""+'"'+str(density)+'"'+"""
	defaultBulkModulus="""+'"'+str(bulk)+'"'+"""
	defaultShearModulus="""+'"'+str(shear)+'"'+"""
	tensileStrength="""+'"'+str(tensileStrength)+'"'+"""
	compressiveStrength="""+'"'+str(compressiveStrength)+'"'+"""
	maximumStrength="""+'"'+str(maximumStrength)+'"'+"""
	crackSpeed="""+'"'+str(crackSpeed)+'"'+"""
	damagedMaterialFrictionSlope="""+'"' + str(damagedMaterialFrictionSlope) + '"' + """
	/>
"""
