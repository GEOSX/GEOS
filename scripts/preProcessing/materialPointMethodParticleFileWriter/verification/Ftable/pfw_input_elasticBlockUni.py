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

sampleWidth = 1.0  # mm
sampleHeight = 1.0 # mm
sampleLength = 1.0 # mm

domainWidth = sampleWidth  # This would be increased for unconfined compression.
domainHeight = sampleHeight
domainLength = sampleLength

pfw["xmin"] = 0.0             # mm
pfw["xmax"] = domainWidth    # mm
pfw["ymin"] =-0.5*domainHeight # mm
pfw["ymax"] = 0.5*domainHeight # mm
pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

refine=2  # partitions in each direction
cpp=3     # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="00:30:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=stopTime
pfw['lastRestartBufferInSeconds'] = 0.

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25  
pfw["initialDt"]=1e-16
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"]=stopTime/2000
pfw["boxAverageHistory"]=1
pfw["boxAverageWriteInterval"]=stopTime/2000

pfw["solverProfiling"]=0         
pfw["frictionCoefficient"]=0.25  

pfw["updateMethod"]="XPIC"
pfw["updateOrder"]=2

# DEFORMATION -----------------------------------------------------------------------------

pfw["fTableInterpType"]='Smoothstep'
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,	 1,	    1,	1],
               [10,	 1.001,	1,	1],
               [20,	 1.002,	1,	1],
               [25,	 1,	    1,	1]
               ]

pfw["prescribedBcTable"]=1    
pfw["bcTable"]=[[0,   2, 2, 2, 2, 2, 2],
                [10,  2, 2, 0, 0, 0, 0],
                [100, 2, 2, 0, 0, 0, 0]]

# MATERIAL PROPERTIES --------------------------------------------------------------------

pfw["materials"] = [ "aluminum" ]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="2.7"
	defaultBulkModulus="70.0"
	defaultShearModulus="24.0"/>
"""

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],v=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]
