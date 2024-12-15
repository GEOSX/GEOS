# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# Test of the prescibed F-table and boundary condition table capability.

# An elastic block is loaded in tensile uniaxial strain, using prescibed F loading.
# The lateral boundaries are then freed using outflow BCs so that subsequent loading is uniaxial stress.

# Material is then unloaded to zero strain.

# Reactions are plotted against the expected result for an elastic material

# Initial slope is constrained modulus, 2nd leg is young's modulus, 3rd leg unloads with young's modulus.

# The Mathematica notebook 'rotateParticles.nb' can be used to rotate particles in-place to see the effect of grid-alignment.

pfw = {}
pfw["runDebug"] = False
stopTime = 25.0

# DOMAIN ---------------------------------------------------------------------------------

cylinderRadius = 2.5    # 2.5 # cylinder radius in mm
endTime = 400                  # 100     		
endStretch = 0.985            # 0.975 

domainWidth = 2.25*cylinderRadius  # This would be increased for unconfined compression.
domainHeight = 2.25*cylinderRadius
domainLength = 2.25*cylinderRadius

pfw["xmin"] = -0.5 * domainWidth             # mm
pfw["xmax"] =  0.5 * domainWidth    # mm
pfw["ymin"] = -0.5 * domainWidth # mm
pfw["ymax"] =  0.5 * domainWidth # mm
pfw["zmin"] =  0.0 # mm
pfw["zmax"] =  domainHeight # mm

refine=2  # partitions in each direction
cpp=10     # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*(cpp+2)  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*(cpp+2)  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="00:03:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/112.)) 
pfw["mSubmitJobs"]=True
pfw["autoRestart"] =True
pfw["mBank"]="__M_BANK__"

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------
pfw["debugFlag"]=0

pfw["endTime"]=endTime
pfw["plotInterval"] = endTime / 40.0
pfw["restartInterval"] = endTime / 12
pfw['lastRestartBufferInSeconds'] = 0.

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.125  
pfw["initialDt"]=1e-16
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"]=stopTime/(2000 * (np.pi/3.141593) )
# pfw["boxAverageHistory"]=0
# pfw["boxAverageWriteInterval"]=stopTime/2000

pfw["solverProfiling"]=0         
pfw["frictionCoefficient"]=0.01  

pfw["updateMethod"]="XPIC"
pfw["updateOrder"]=2
pfw["useAPIC"]=0
pfw["useInteralForceAsFaceReaction"]=1

# DEFORMATION -----------------------------------------------------------------------------

pfw["fTableInterpType"]='Smoothstep'
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,	 1,	    1,	1],
               # [10,	 1.0125,	1,	1],
               # [20,	 1.025,	1,	1],
               [0.5 * endTime,	 1,	    1,	endStretch],
               [1.0 * endTime,	 1,	    1,	1]
               ]

pfw["boundaryConditionTypes"] = "{ 2, 2, 2, 2, 2, 2 }"
# pfw["prescribedBcTable"]=1    
# pfw["bcTable"]=[[0,   2, 2, 2, 2, 2, 2],
#                 # [10,  2, 2, 0, 0, 0, 0],
#                 [ 30 ,  2, 2, 0, 0, 0, 0],
#                 [100, 2, 2, 0, 0, 0, 0]]

# MATERIAL PROPERTIES --------------------------------------------------------------------

pfw["materials"] = [ "aluminum" ]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="1.935"
	defaultBulkModulus="0.13889"
	defaultShearModulus="0.10417"/>
"""

# GEOMETRY OBJECTS -------------------------------------------------------

# block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],
     # v=[0.0,0.0,0.0],mat=0,group=0)
block = geom.cylinder('cylinder',[0.0, 0.0, 0.0], [0.0, 0.0, pfw["zmax"]], 
	cylinderRadius, v=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# #------------------------------------------------------------------------------#
# thisCPP = 12  # cell per partition count
# cylinderRadius = 2.5    # 2.5 # cylinder radius in mm
# endTime = 25                  # 100     		
# endStretch = 0.99            # 0.975 

# # sample material parameters
# sampleDensity = 1.935         # 1.935 #mg/mm3
# sampleYoungsModulus = 0.25         # 0.25 # GPa
# samplePoissonRatio = 0.2         # 0.2   # 0.4

# thisRefine = 2  # domain partitioning 
# sampleBulk = sampleYoungsModulus / (3*(1-2*samplePoissonRatio))
# sampleShear = sampleYoungsModulus / (2*(1+samplePoissonRatio))

# pfw["cflFactor"]=0.25  

# # Domain ----------------------------------------------------------------------#
# pfw["xpar"] = 1*(thisRefine)  	# grid partitions
# pfw["ypar"] = 1*(thisRefine)
# pfw["zpar"] = 1*(thisRefine)

# cppZ=thisCPP
# cppXY = cppZ+0#1
# pfw["nI"]=pfw["xpar"]*cppXY  # grid cells in the x-direction
# pfw["nJ"]=pfw["ypar"]*cppXY  # grid cells in the y-direction
# pfw["nK"]=pfw["zpar"]*cppZ 	# grid cells in the z-direction
# pfw["ppc"]=2   						# particles per cell in each direction

# # BATCH PARAMETERS  --------------------------------------------------------

# pfw["mBatch"]=True
# pfw["mWallTime"]="00:30:00"
# pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
# pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/112.)) 
# pfw["mSubmitJobs"]=False

# # Define all the geometric objects --------------------------------------------#
# domainHeight = 1.0 * (2.0 * cylinderRadius)	# 1:1:1 ratio domain
# domainWidth = pfw["nI"]/pfw["nK"]*domainHeight

# pfw["xmin"] = -0.5*domainWidth		# mm
# pfw["xmax"] =  0.5*domainWidth		# mm
# pfw["ymin"] = -0.5*domainWidth		# mm
# pfw["ymax"] =  0.5*domainWidth		# mm
# pfw["zmin"] = 0.0 					# mm
# pfw["zmax"] = domainHeight  		# mm

# dx = (pfw["xmax"]-pfw["xmin"])/(pfw["nI"]-2)/pfw["ppc"]
# dy = (pfw["ymax"]-pfw["ymin"])/(pfw["nJ"]-2)/pfw["ppc"]
# dz = (pfw["zmax"]-pfw["zmin"])/(pfw["nK"]-2)/pfw["ppc"]

# # GEOSX MPM input parameters --------------------------------------------------#

# pfw["objects"] = []

# sample = geom.cylinder('cylinder',[0.0, 0.0, 0.0], [0.0, 0.0, domainHeight], 
# 	cylinderRadius, v=[0.0,0.0,0.0],mat=0,group=0)

# # sample = geom.box('cylinder',[0.0-cylinderRadius, 0.0-cylinderRadius, 0.0], 
# # 	[0.0+cylinderRadius, 0.0+cylinderRadius, domainHeight], 
# # 	v=[0.0,0.0,0.0],mat=0,group=0)
# # basePlaten = geom.box('base', [-100, -100, -100], [100, 100, basePlatenDepth], 
# # 	v=[0,0,0], mat=0, group=0)

# pfw["objects"] = pfw["objects"] + [sample]

# # Solver Values ----------------------------------------------------------------


# # Process values:
# #------------------------------------------------------------------------------#
# pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
# pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/112.)) 
# #------------------------------------------------------------------------------#

# # Solver Values ----------------------------------------------------------------
# pfw["endTime"] = endTime
# pfw["plotInterval"] = endTime / 12.0
# pfw["restartInterval"] = endTime / 1
# pfw["reactionHistory"] = 1
# pfw["reactionWriteInterval"] = endTime / 500.0
# pfw["boxAverageHistory"] = 0

# pfw["updateMethod"]="XPIC"
# pfw["updateOrder"]=2

# pfw["timeIntegrationOption"]="ExplicitDynamic"
# pfw["cflFactor"]=0.5   
# pfw["initialDt"]=1e-16

# # pfw["cpdiDomainScaling"] = 1
# # pfw["bodyForce"] = [ 0, 0, 0 ]
# # pfw["boundaryConditionTypes"] = "{ 2, 2, 2, 2, 2, 2 }"

# pfw["prescribedBcTable"]=1    
# pfw["bcTable"]=[[0,   2, 2, 2, 2, 2, 2],
#                 # [10,  2, 2, 0, 0, 0, 0],
#                 [20,  2, 2, 0, 0, 0, 0],
#                 [100, 2, 2, 0, 0, 0, 0]]

# pfw["fTableInterpType"]='Smoothstep'
# pfw["prescribedBoundaryFTable"]=1
# pfw["fTable"]=[[0,	 1,	    1,	1],
#                # [10,	 1.0125,	1,	1],
#                # [20,	 1.025,	1,	1],
#                [25,	 1,	    1,	0.99]
#                ]

# # pfw["damageFieldPartitioning"] = 0
# # pfw["needsNeighborList"] = 0
# pfw["useDamageAsSurfaceFlag"] = 0
# pfw["frictionCoefficient"] = 0.25

# pfw["materials"] = [ "LESample" ]
# pfw["materialPropertyString"]="""
# <ElasticIsotropic
# 	name="LESample"
# 	defaultDensity="""+'"'+str(sampleDensity)+'"'+"""
# 	defaultBulkModulus="""+'"'+str(sampleBulk)+'"'+"""
# 	defaultShearModulus="""+'"'+str(sampleShear)+'"'+"""/>
# """

# Define pfw["particleFileFields"]
# #------------------------------------------------------------------------------#
# pfw["particleFileFields"] = ["Velocity", "MaterialType", "ContactGroup", 
# 							 "SurfaceFlag", "StrengthScale", "RVector"]

#------------------------------------------------------------------------------#
# pfw["plottableFields"] =["particleReferencePosition", "particleDensity",
#  "particleMaterialType", "particleCenter", "particleVelocity", "particleVolume",
#  "particleID", "particleAcceleration", "particleStress", "particleDamage"]
