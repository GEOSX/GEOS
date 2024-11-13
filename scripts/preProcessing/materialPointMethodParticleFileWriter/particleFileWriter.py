# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 09:00:00 2017

@author: homel1
         crook5

Script to automate generation of a GEOSX-MPM input file and particle file
based on geometric objects.

Particle file format depends on specified plotting fields.

"""
from __future__ import print_function    # (at top of module)
from __future__ import division
from __future__ import unicode_literals
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree # nearest neighbor search with KDTree
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import datetime
import time
import datetime                               # used for date stamp
import os                                     # operating sys commands, pwd, etc.
from subprocess import call                   # lets you call shell commands, i.e. call(["ln", "-s",".","run_dir"])
import subprocess                             # lets you call msub and get jobid
import sys                                    # to access command arguments.
import importlib
import random                                 # used to define a random material type
import pfw_geometryObjects as geom            # this contains all the geometry object functions for pfw
import math
import getpass
import platform
import difflib

# ============================================================================================
# MACHINE-SPECIFIC Calculations.
# ============================================================================================

# List of cores per node of LC (or other) machines.  TODO: allow user to append this in user defs.
machineList = {
  'lassen':44,
  'dane':112,
  'ruby':56,
  'rzhound':56,
  'tioga':64
}

node = platform.node()
for key, value in machineList.items():
  if key in node:
    machine = key
    coresPerNode = value
    # This could be unsafe if someone added a machines name we use elsewhere.
    # Currently we test if lassen=True for various MPI tasks.
    exec(key+'=True')
  else:
    exec(key+'=False')

# # MPI specific variables
# there seems to be an issue with mpi4py and subprocess launching 
# slurm scripts that include srun. We've included the '#SBATCH --export=NONE' command
# in launched scripts below, which seems to fix it, but be aware there may be
# issues when launching subprocesses from this script, even if no MPI commands are used
# or if num_ranks=0.
#
if lassen:
  rank = 0
  num_ranks = 1
else:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD # gets communication pool 
  rank = comm.Get_rank()  # gets rank of current process 
  num_ranks = comm.Get_size() # total number of processes


# ============================================================================================
# BEGIN FUNCTION DEFINITIONS
# ============================================================================================

#This code calculates the similarity between two strings using the ndiff method from the difflib library. 
def compute_similarity(input_string, reference_string):
#The ndiff method returns a list of strings representing the differences between the two input strings.
    diff = difflib.ndiff(input_string, reference_string)
    diff_count = 0
    for line in diff:
      # a "-", indicating that it is a deleted character from the input string.
        if line.startswith("-"):
            diff_count += 1
# calculates the similarity by subtracting the ratio of the number of deleted characters to the length of the input string from 1
    return 1 - (diff_count / len(input_string))

# ============================================================================================
# END FUNCTION DEFINITIONS
# ============================================================================================

# There are places where we compare numpy array to string, currently works, might not in the future
# fix it, but for now so we don't spam the log file:
# https://stackoverflow.com/questions/40659212/futurewarning-elementwise-comparison-failed-returning-scalar-but-in-the-futur
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# ===========================================
# READ FROM USER INPUT FILE
# ===========================================

username = getpass.getuser()
userDefsFile = str('userDefs_'+str(username))

try:
  f = open(userDefsFile+".py")
  # Do something with the file
except IOError:
  print("Please create a version of userDefs_OUN.py consistent with your user name and geosx location")
finally:
  f.close()

userDefs = importlib.import_module(userDefsFile)
geosPath = userDefs.geosPath

# read timeStamp as argument.  If it's missing give a warning.
if len(sys.argv) < 2:
  print("usage : particleFileWriter.py fileName \nYou must specify the fileName as arguments")
  sys.exit(1)

# inputFile name (strip the extension just in case it was called that way)
inputFile = (str(sys.argv[1])).replace('.py',"")

# Temporarily check if you want to regenerate the particle file
generateParticleFile = True
if len(sys.argv) > 2:
  generateParticleFile = bool(int(sys.argv[2]))
print("Generate particle file? ", generateParticleFile)

#import inputFile as job    # Input parameters in separate file.
job = importlib.import_module(inputFile)

# file name prefix is date stamped
now = datetime.datetime.now()
timeStamp = str(now.year)[-2:]+str(now.month).zfill(2)+str(now.day).zfill(2) + \
            str(now.hour).zfill(2)+str(now.minute).zfill(2)+str(now.second).zfill(2)

# Current working directory
PWD = os.getcwd()

# Variable Info
# sortObjects: This is useful for large simulations with many object.  It constructs a subset of objects
#              for each slice, and only searches those. For granular asseblies of objects this should be
#              much faster.  Not all objects have the necessary xmin and xmax attribute, so this isnn't
#              default behavior.  This shouldn't change order of objects for first-in priority.


# If a specific variable needs a check before being written to solver string a handle to that function is added as a value in the dictionary below
# Value contains ( default value, flag to include in xml mpm solver parameter string if not specified or not )
parameters = { 'runDebug' : ( False, False ),
               'stopTime' : ( False, False),
               'mBank' : ( None, False ),
               'mWallTime' : ( "00:30:00", False ),
               'mBatch' : ( True, False ),
               'mCores' : ( 1, False ),
               'mNodes' : ( 1, False ), 
               'mSubmitJobs' : ( False, False ),
               'autoRestart' : ( False, False ),
               'mPartition' : ( "pbatch", False ),
               'periodic' : ( [False, False, False], False ),
               'xpar' : ( 1, False ),
               'ypar' : ( 1, False ),
               'zpar' : ( 1, False ),
               'nI' : ( 5, False ),
               'nJ' : ( 5, False ),
               'nK' : ( 5, False ),
               'ppc' : ( 2, False ),
               'ppcx' : ( None, False ),
               'ppcy' : ( None, False ),
               'ppcz' : ( None, False ),
               'xmin' : ( -0.5, False ),
               'xmax' : ( 0.5, False ),
               'ymin' : ( 0.0, False ),
               'ymax' : ( 1.0, False ),
               'zmin' : ( -0.5, False ),
               'zmax' : ( 0.5, False ),
               'lastRestartBufferInSeconds' : ( 10, False ),
               'objects' : ( None, False ),
               'sortObjects' : ( False, False ),
               'timeIntegrationOption': ( None, True ),
               'updateMethod': ( None, True ),
               'updateOrder': ( None, True ),
               'cflFactor': ( None, True ),
               'initialDt': ( None, True ),
               'solverProfiling': ( None, True ),
               'cpdiDomainScaling': ( None, True ),
               'damageFieldPartitioning': ( None, True ), 
               'planeStrain': ( False, True ),
               'needsNeighborList': ( None, True ),
               'reactionHistory': ( None, True ),
               'reactionWriteInterval': ( None, True ),
               'boxAverageHistory': ( None, True ),
               'boxAverageWriteInterval': ( None, True ),
               'prescribedBcTable': ( None, True ),
               'prescribedFTable': ( None, True ),
               'prescribedBoundaryFTable': ( None, True ),
               'fTable': ( None, True ),
               'fTableInterpType': ( None, True ),
               'stressControl': ( None, True ),
               'stressTable': ( None, True ),
               'stressControlKp': ( None, True ),
               'stressControlKi': ( None, True ),
               'stressControlKd': ( None, True ),
               'stressTableInterpType': ( None, True ),
               'boundaryConditionTypes': ( None, True ),
               'bcTable': ( None, True ),
               'useEvents': ( None, True ),
               'bodyForce' : ( None, True ),
               'generalizedVortexMMS' : ( None, True ),
               'debugFlag' : ( None, True ),
               'frictionCoefficient' : ( None, True ),
               'frictionCoefficientTable' : ( None, True ),
               'frictionCoefficientRuleOfMixtures' : ( None, True ),
               'useDamageAsSurfaceFlag' : ( False, True ),
               'neighborRadius' : ( None, True ),
               'minParticleJacobian' : ( None, True ),
               'maxParticleJacobian' : ( None, True ), 
               'debugFlag' : ( None, True ),
               'materials' : ( None, False ),
               'materialPropertyString' : ( None, False ),
               'endTime' : ( 1.0, False ),
               'plotInterval' : ( None, False ),
               'restartInterval' : ( None, False ),
               'useSinusoidalDamageField' : ( False, False),
               'wavyCrack' : ( False, False),
               'particleRefinement' : ( None, False ),
               'mpmEventsString' : ( '', False ),
               'maxParticleVelocity' : ( 10.0, True ),
               'cohesiveFieldPartitioning' : ( 0, True),
               'enableCohesiveFailure' : ( 0, True ),
               'maxCohesiveNormalStress' : ( 0.01, True ),
               'maxCohesiveShearStress' : ( 0.01, True ),
               'characteristicNormalDisplacement' : ( 0.01, True ),
               'characteristicTangentialDisplacement' : ( 0.01, True ),
               'maxCohesiveNormalDisplacement' : ( 0.01, True ),
               'maxCohesiveTangentialDisplacement' : ( 0.01, True ),
               'prescribedBoundaryTransverseVelocities' : ( [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]], False ),
               'particleFileFields' : (["Velocity",  # +3
                                        "MaterialType", # +1
                                        "ContactGroup", # +1
                                        "SurfaceFlag", # +1
                                        "RVector"] , False) } # +9

# This must correspond to the order in which fields are written in the particle file below
particleFieldOrder=["Velocity",
                    "MaterialType",
                    "ContactGroup",
                    "SurfaceFlag",
                    "Damage",
                    "Porosity",
                    "Temperature",
                    "StrengthScale",
                    "RVector",
                    "MaterialDirection",
                    "SurfaceNormal",
                    "SurfacePosition",
                    "SurfaceTraction",
                    "ShrinkageFlag"]


# New dictionary input approach
pfw = job.pfw if hasattr(job, 'pfw') else {}



tabIndent = 3*"  "  
parameterStrings = []
for paramName, paramTuple in parameters.items():
  # Check if param is passed from input script, if not assign default value
  paramValue = paramTuple[0]
  if paramName in pfw:
    paramValue = pfw[paramName] 
    if paramValue != None and paramTuple[1]:
      parameterStrings.append(tabIndent + paramName + '="' + str(paramValue).replace('[','{').replace(']','}') + '"' + '\n')

  # Add global variable to be used by particle file writer
  globals()[paramName] = paramValue 

# Add all remaining variables to mpmSolverParameterString that aren't specified here, but don't add them as global variables for script
for paramName, paramValue in pfw.items():
  if paramName not in parameters:

    maxDiff = 0.0
    closestParam = ""
    for p, paramTuple in parameters.items():
      newDiff = compute_similarity( paramName, p )
      if newDiff > 0.5 and newDiff > maxDiff:
        maxDiff = newDiff
        closestParam = p

    if maxDiff != 0.0:
      print( paramName + " parameter not found, did you mean " + closestParam +  "?")
    else:
      print( paramName + " parameter not found")

    parameterStrings.append(tabIndent + paramName + '="' + str(paramValue).replace('[','{').replace(']','}').replace('\'','') + '"' + '\n')
    

# Remove new line from last parameter to be added (for xml formatting)
parameterStrings[-1] = parameterStrings[-1].replace('\n','')
mpmSolverParameterString = ''.join(parameterStrings)

# Remove fields form particleFileOrder the user does not wish to write to particle field (particleFieldOrder preseveres the correct header order for values written to the particle file)
particleFieldOrder = [ f for f in particleFieldOrder if f in particleFileFields ]

# Interior Domain
# The geosx input xml file specifies and xmin,xmax, ni, etc. based on the total domain
# size (with ghost cells).  The user inputs in the python input file will be physical
# domain size (interior only) and actual number of grid cells, to make partitioning
# easier.  here we calculate the true domain extents and particle sizes.
# Discretization
NI = nI   # total grid cells in the x-direction (counting ghosts)
NJ = nJ   # total grid cells in the y-direction (counting ghosts)
NK = nK   # total grid cells in the z-direction (counting ghosts)
ppcx = ppcx if ppcx != None else ppc
ppcy = ppcy if ppcy != None else ppc
ppcz = ppcz if ppcz != None else ppc

# list of geometry creation objects
if objects == None:
  objects = job.make_objects()

# Batch parameters for GEOSX runs.  An error will result if there are too many cores for
# a low resolution simulation.  If there is insufficient run-time to obtain a signal
# for a given run, that run will have its results ommited from the Hugoniot analysis.
mPartition = mPartition if not runDebug else "pdebug" 

# We used to set this manually, but coresPerNode changes with each machine.
# This will ensure consistency since we now have that value for each platform.
mNodes= int(np.ceil(float(mCores)/float(coresPerNode))) 
print('machine = ',machine,', mNodes = ',mNodes,', mCores = ',mCores,', coresPerNode = ',coresPerNode)


if mBank == None:
  # Get default bank from userdefs
  username = getpass.getuser()
  userDefsFile = str('userDefs_'+str(username))
  userDefs = importlib.import_module(userDefsFile)
  mBank = userDefs.defaultBank




[wH,wM,wS]=mWallTime.split(":")
wallTimeMinutes=int(wH)*60+int(wM)
if wallTimeMinutes > 60 and runDebug:
  print("Wall time of debug job exceeded 60 minutes and was reset")
  wallTimeMinutes = 60
  mWallTime="01:00:00"

maxRestartTime = wallTimeMinutes*60 - min(wallTimeMinutes,lastRestartBufferInSeconds)
mWallTimeMinutes=str(wallTimeMinutes)

coreHours = float(int(mCores))*(float(int(wH))+float(int(wM))/60.+float(int(wS))/3600.)
cpuTimeCost = (23295./2000000.)*coreHours
print('Approximate LC bank time cost for this simulation is $',cpuTimeCost,'.')

# Material array
if materials == None:
  # Throw error
  print( "A materials array must be specified in the input file!" )
  comm.abort()

matsOrig = materials
mats = str(matsOrig).replace("[",'"{')
mats = mats.replace("]",'}"')
mats = mats.replace("'","")

particleTypesPerMat = [set() for m in materials]

particleRefinement = particleRefinement if particleRefinement != None else [ 1 for i in range(len(materials)) ]# Create list of size materials all ones

# ============================================================================================
# ERROR Checking
# ============================================================================================
if rank == 0:
  no_errors = True

  if ( planeStrain and NK != 3):
    print('nK should = 3 for plane Strain!!')
    no_errors = False

  if ( planeStrain and zpar != 1):
    print('zpar should = 1 for plane Strain!!')
    no_errors = False

  if not no_errors:
    if not lassen:
      comm.Abort()
    else:
      sys.exit()

# ============================================================================================
# CREATE PARTICLE FILE
# ============================================================================================

timer = time.time()
particleFileName = 'mpmParticleFile_'+inputFile.replace('pfw_input_',"")

# interior discretizatiom
nI = NI if periodic[0] else NI-2   # interior grid cells in the x-direction
nJ = NJ if periodic[1] else NJ-2   # interior grid cells in the y-direction
nK = NK if periodic[2] else NK-2   # interior grid cells in the z-direction

# grid cell spacing
dX = (xmax - xmin)/nI
dY = (ymax - ymin)/nJ
dZ = (zmax - zmin)/nK

# total domain dimensions (counting ghosts)
XMIN = xmin if periodic[0] else xmin - dX 
XMAX = xmax if periodic[0] else xmax + dX
YMIN = ymin if periodic[1] else ymin - dY
YMAX = ymax if periodic[1] else ymax + dY
ZMIN = zmin if periodic[2] else zmin - dZ
ZMAX = zmax if periodic[2] else zmax + dZ

# particles in each direction
ni = ppcx*nI
nj = ppcy*nJ
nk = 1 if planeStrain else ppcz*nK

# particle dimensions
dx = dX/ppcx
dy = dY/ppcy
dz = dZ if planeStrain else dZ/ppcz


# Sum particle volume to compute volume fraction.
particleVolume = 0.
domainVolume = ( pfw["xmax"] - pfw["xmin"] )*( pfw["ymax"] - pfw["ymin"] )*( pfw["zmax"] - pfw["zmin"] ) 

if generateParticleFile:
  print('Writing particle file...')

  # Delimiter for particle file
  delim = '\t'

  rank_particleFileName = particleFileName + '_' + str(rank)
  particleFile = open(rank_particleFileName, 'w')

  if(planeStrain):
    surfaceDepth = 2*np.sqrt(dX*dX + dY*dY) / min(ppcx, ppcy)
  else:
    surfaceDepth = np.sqrt(dX*dX + dY*dY + dZ*dZ) / min(ppcx, min( ppcy, ppcz))

  # list of i,j,k indices at the center of each partition.
  ipc=np.empty(xpar,dtype=int)
  jpc=np.empty(ypar,dtype=int)
  kpc=np.empty(zpar,dtype=int)
  for p in range(xpar):
    ipc[p] = np.ceil( (p+0.5)*ppcx*NI/xpar ) - ppcx
  for p in range(ypar):
    jpc[p] = np.ceil( (p+0.5)*ppcy*NJ/ypar ) - ppcy
  for p in range(zpar):
    if planeStrain:
      kpc[p] = 1
    else:
      kpc[p] = np.ceil( (p+0.5)*ppcz*NK/zpar ) - ppcz

  if rank == 0:
    print('xpar = ',xpar,', ppcx*NI = ',ppcx*NI,', partition centers at i = ',ipc)
    print('ypar = ',ypar,', ppcy*NJ = ',ppcy*NJ,', partition centers at j = ',jpc)
    print('zpar = ',zpar,', ppcz*NK = ',ppcz*NK,', partition centers at k = ',kpc)

    # loop through domain and fill particles.

    print("ni=",ni,"nj=",nj,",nk=",nk)

  # Make sure no errors are detected on rank 0 before proceeding with particle generation
  if not lassen:
    comm.Barrier()

  max_num_particles = ni*nj*nk

  n_p = 0
  for i in range(1+math.ceil(ni*rank/num_ranks), math.ceil(ni*(rank+1)/num_ranks)+1):
    print('creating files, row',i,'/',ni,', estimated time remaining = ',(ni-i)/i*(time.time()-timer),'s')
    x = xmin + (i - 0.5)*dx

    if(sortObjects):
      sliceObjects = [ob for ob in objects if (ob.xMin() <= x and ob.xMax() >= x)]
    else:
      sliceObjects = objects

    for j in range(1,nj+1):
        y = ymin + (j - 0.5)*dy
        for k in range(1,nk+1):
          z = zmin + (k - 0.5)*dz

          match = False
          ob = 0

          while (ob < len(sliceObjects) and match == False):
            object = sliceObjects[ob]
            ob = ob+1
            pt = np.array([x,y,z])
            voxelFlag = object.isInterior( pt, surfaceDepth ) # Voxel flags greater than or equal to 0 denote interior regions of object
            if( voxelFlag >= 0 ):
              match = True

              # Get the particle type, default is CPDI
              if hasattr(object, 'getParticleType'):
                particleType = object.getParticleType( pt )
              else:
                particleType = object.particleType( object, pt ) if callable( object.particleType ) else object.particleType

              # set the particle damage to be the object damage unless
              # overwritten by surface flags.
              if "Velocity" in particleFileFields:
                if hasattr(object, 'getVelocity'):
                  velocity = object.getVelocity( pt )
                else:
                  velocity = object.v( object, pt ) if callable( object.v ) else object.v

              if "Damage" in particleFileFields:
                # set the particle damage to be the object damage unless
                # overwritten by surface flags.
                if hasattr(object, 'getDamage' ): # some internal function defines spatially varing initial damage
                  damage = object.getDamage( pt )
                else: # damage is constant for the object:
                  damage = object.damage(object, pt ) if callable( object.damage ) else object.damage

                # force damage to be a scalar value.
                if(useSinusoidalDamageField):
                  damage = 0.5*( np.sin( 2.0*np.pi*(x-xmin)/(xmax-xmin) )*np.sin(  3.0*np.pi*(y-ymin)/(ymax-ymin) ) + 1 )

                if(wavyCrack):
                  h = 2*dy
                  val = 0.25 * (ymax - ymin) * np.cos(2.0*3.1415926535*(x-xmin)/(xmax-xmin)) + ymin + 0.5 * (ymax - ymin)
                  if( val - h < y and y < val + h ):
                    damage = 1
                  else:
                    damage = 0

              if "Porosity" in particleFileFields:
                if hasattr(object, 'getPorosity' ):
                  porosity = object.getPorosity( pt )
                else: # damage is constant for the object:
                  porosity = object.porosity(object, pt ) if callable( object.porosity ) else object.porosity

              if "Temperature" in particleFileFields:
                if hasattr(object, 'getTemperature' ):
                  temperature = object.getTemperature( pt )
                else: # temperature is constant for the object:
                  temperature = object.temperature(object, pt ) if callable( object.temperature ) else object.temperature


              # Particle Surface Flags:
              #  enum struct SurfaceFlag : integer
              #  {
              #    0: Interior,
              #    1: FullyDamaged,
              #    2: Surface,
              #    3: Cohesive    
              #  };


              if "SurfaceFlag" in particleFileFields:
                # Surface flags need to be 2 for polymer model to not crash
                surfaceFlag = voxelFlag #object.isSurface( pt, surfaceDepth ) 

              surfacePosition = np.array([0.0, 0.0, 0.0])
              if "SurfacePosition" in particleFileFields and surfaceFlag != 0:
                 surfacePosition = object.getSurfacePosition( pt )
                
              surfaceNormal = np.array([0.0, 0.0, 0.0])
              if "SurfaceNormal" in particleFileFields and surfaceFlag != 0:
                surfaceNormal = object.getSurfaceNormal( pt )

              surfaceTraction = np.array([0.0,0.0,0.0])  
              if "SurfaceTraction" in particleFileFields and surfaceFlag != 0:
                if hasattr(object, 'getSurfaceTraction'):
                  surfaceTraction = object.getSurfaceTraction( pt )
                else:
                  surfaceTraction = object.surfaceTraction(object, pt) if callable(object.surfaceTraction) else object.surfaceTraction

              if "MaterialType" in particleFileFields:
                if hasattr(object, 'getMat'):
                  mat = int( object.getMat( pt ) )
                else:
                  mat = int( object.mat(object, pt) if callable(object.mat) else object.mat )

              if "ContactGroup" in particleFileFields:
                if hasattr(object,'getGroup'):
                    group = object.getGroup( pt )
                else:
                    if hasattr(object, 'group'):
                      group = object.group
                    else:
                      group = 0

              if "MaterialDirection" in particleFileFields:
                # material direction, These will be read from object, will default to [1,0,0] if not 
                # specified.
                if hasattr( object, 'getMatDir' ):
                  matDir = object.getMatDir( pt )
                else:
                  if hasattr( object, 'matDir' ):
                    matDir = object.matDir( object, matDir ) if callable( object.matDir ) else object.matDir
                  else:
                    matDir = np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]])
  
              if "StrengthScale" in particleFileFields:
                # object may have some constant or spatially varying strength scale.
                # which is used by some material properties to modify strength.
                # this is an alternative to using internal weibull scaling, but may
                # be useful for defining other spatial distributions or object
                # specific weak layers, etc. defaults to 1 if not specified.   
                if hasattr(object,'getStrengthScale'):
                    strengthScale = object.getStrengthScale( pt )
                else:
                    if hasattr(object, 'strengthScale'):
                      strengthScale = object.strengthScale( object, pt ) if callable(object.strengthScale) else object.strengthScale
                    else:
                      strengthScale = 1.0

              # Check for fields with None values
              # Get the particle type, default is CPDI
              if particleType == None:
                print("ParticleType value from", object.name,"was None!")
                sys.exit(0)

              if "Velocity" in particleFileFields and all(v is None for v in velocity):
                print("Velocity value from", object.name,"was None!")
                sys.exit(0)

              if "Damage" in particleFileFields and damage == None:
                print("Damage value from", object.name,"was None!")
                sys.exit(0)

              if "Porosity" in particleFileFields and porosity == None:
                print("Porosity value from", object.name,"was None!")
                sys.exit(0)

              if "Temperature" in particleFileFields and temperature == None:
                print("Temperature value from", object.name,"was None!")
                sys.exit(0)

              if "SurfaceFlag" in particleFileFields and surfaceFlag == None:
                print("SurfaceFlag value from", object.name,"was None!")
                sys.exit(0)

              if "SurfacePosition" in particleFileFields and surfaceFlag != 0 and all(v is None for v in surfacePosition):
                print("SurfacePosition value from", object.name,"was None!")
                sys.exit(0)

              if "SurfaceNormal" in particleFileFields and surfaceFlag != 0 and all(v is None for v in surfaceNormal):
                print("SurfaceNormal value from", object.name,"was None!")
                sys.exit(0)

              if "SurfaceTraction" in particleFileFields and surfaceFlag != 0 and all(v is None for v in surfaceTraction):
                print("SurfaceTraction value from", object.name,"was None!")
                sys.exit(0)

              if "MaterialType" in particleFileFields and mat == None:
                print("MaterialType value from", object.name,"was None!")
                sys.exit(0)

              if "ContactGroup" in particleFileFields and group == None:
                print("ContactGroup value from", object.name,"was None!")
                sys.exit(0)

              if "MaterialDirection" in particleFileFields and all(v is None for v in matDir):
                print("MaterialDirection value from", object.name,"was None!")
                sys.exit(0)
  
              if "StrengthScale" in particleFileFields and strengthScale == None:
                print("Strengthscale value from", object.name,"was None!")
                sys.exit(0)

              particleTypesPerMat[mat].add( particleType )

              dxr = dx/particleRefinement[mat]
              dyr = dy/particleRefinement[mat]
              dzr = dz if planeStrain else dz/particleRefinement[mat]
              
              for ri in range(particleRefinement[mat]):
                for rj in range(particleRefinement[mat]):
                  for rk in range( ( 1 if planeStrain else particleRefinement[mat] ) ):
                    xr = x - 0.5*dx + (ri + 0.5)*dxr
                    yr = y - 0.5*dy + (rj + 0.5)*dyr
                    zr = z - 0.5*dz + (rk + 0.5)*dzr
                    
                    # only create particle if mat>=0, that way we can specify mat=-1 to generate a 
                    # defect or void.
                    if(mat>=0):
                      n_p = n_p + 1

                      particleVolume += dxr*dyr*dzr

                      pString = str(n_p) + delim \
                                + str(xr) + delim \
                                + str(yr) + delim \
                                + str(zr) + delim \
                                + str(particleType)

                      if "Velocity" in particleFileFields:
                        pString = pString + delim + str(velocity[0]) + delim \
                                                  + str(velocity[1]) + delim \
                                                  + str(velocity[2])

                      if "MaterialType" in particleFileFields:
                        pString = pString + delim + str(mat)
                        
                      if "ContactGroup" in particleFileFields:
                        pString = pString + delim + str(group)

                      if "SurfaceFlag" in particleFileFields:
                        pString = pString + delim + str(surfaceFlag)

                      if "Damage" in particleFileFields:
                        pString = pString + delim + str(damage)

                      if "Porosity" in particleFileFields:
                        pString = pString + delim + str(porosity)

                      if "Temperature" in particleFileFields:
                        pString = pString + delim + str(temperature)

                      if "StrengthScale" in particleFileFields:
                        pString = pString + delim + str(strengthScale)

                      if "RVector" in particleFileFields:
                        pString = pString + delim + str(dxr*0.5) + delim \
                                                  + str(0) + delim \
                                                  + str(0) + delim \
                                                  + str(0) + delim \
                                                  + str(dyr*0.5) + delim \
                                                  + str(0) + delim \
                                                  + str(0) + delim \
                                                  + str(0) + delim \
                                                  + str(dzr*0.5)

                      if "MaterialDirection" in particleFileFields:
                        pString = pString + delim + str(matDir[0][0]) + delim \
                                                  + str(matDir[0][1]) + delim \
                                                  + str(matDir[0][2]) 

                      # if "MaterialDirection" in particleFileFields:
                      #   pString = pString + delim + str(matDir[0][0]) + delim \
                      #                             + str(matDir[0][1]) + delim \
                      #                             + str(matDir[0][2]) + delim \
                      #                             + str(matDir[1][0]) + delim \
                      #                             + str(matDir[1][1]) + delim \
                      #                             + str(matDir[1][2]) + delim \
                      #                             + str(matDir[2][0]) + delim \
                      #                             + str(matDir[2][1]) + delim \
                      #                             + str(matDir[2][2])

                      if "SurfaceNormal" in particleFileFields:
                        pString = pString + delim + str(surfaceNormal[0]) + delim \
                                                  + str(surfaceNormal[1]) + delim \
                                                  + str(surfaceNormal[2])
                      
                      if "SurfacePosition" in particleFileFields:
                        pString = pString + delim + str(surfacePosition[0]) + delim \
                                                  + str(surfacePosition[1]) + delim \
                                                  + str(surfacePosition[2])

                      if "SurfaceTraction" in particleFileFields:
                        pString = pString + delim + str(surfaceTraction[0]) + delim \
                                                  + str(surfaceTraction[1]) + delim \
                                                  + str(surfaceTraction[2])

                      if "ShrinkageFlag" in particleFileFields:
                        pString = pString + delim + str( object.flag if hasattr( object, 'flag') else 0)

                      pString = pString +'\n'
                      particleFile.write(pString)

  particleFile.close()

  # Wait for all mpi processes to finish generating particles
  if lassen:
    print('Lassen does not mpi')
  else:
    comm.Barrier()

  n_p = 0
  if rank == 0:
    # Merge each process's particle files together
    # TODO: may be faster to use linux command
    with open(particleFileName, 'w') as outfile:
      # Write column headers
      columnNames = "ID" + delim + "PositionX" + delim + "PositionY" + delim + "PositionZ" + delim + "ParticleType"
      for field in particleFieldOrder:
        fieldString = field

        if field == "Velocity" or field == "SurfaceNormal" or field == "SurfacePosition" or field == "SurfaceTraction" or field == "MaterialDirection":
          fieldString = field + "X" + delim + field + "Y" + delim + field + "Z"

        if field == "RVector":
          fieldString =  field + "XX" + delim + \
                         field + "XY" + delim + \
                         field + "XZ" + delim + \
                         field + "YX" + delim + \
                         field + "YY" + delim + \
                         field + "YZ" + delim + \
                         field + "ZX" + delim + \
                         field + "ZY" + delim + \
                         field + "ZZ"

        # if field == "Velocity" or field == "SurfaceNormal" or field == "SurfacePosition":
        #   fieldString = field + "X" + delim + field + "Y" + delim + field + "Z"

        # if field == "RVector" or field == "MaterialDirection":
        #   fieldString =  field + "XX" + delim + \
        #                  field + "XY" + delim + \
        #                  field + "XZ" + delim + \
        #                  field + "YX" + delim + \
        #                  field + "YY" + delim + \
        #                  field + "YZ" + delim + \
        #                  field + "ZX" + delim + \
        #                  field + "ZY" + delim + \
        #                  field + "ZZ"

        columnNames = columnNames + delim + fieldString

      columnNames = columnNames + "\n"
      outfile.write(columnNames)

      for r in range(0,num_ranks):
        with open(particleFileName + '_'  +str(r)) as infile:
          for line in infile:
            n_p = n_p + 1
            # replace the per-rank id with an ordered global particle id:
            outfile.write(str(n_p) + ' ' + line.split(delim, 1)[1])
        call(["rm",particleFileName+'_'+str(r)],cwd=PWD)

    print('Created n_p = ',n_p,' particles in t = ',time.time()-timer)
    print('Particle volume = ',particleVolume,', Domain volume = ',domainVolume )
    print('Particle volume fraction = ',particleVolume/domainVolume )


# ===========================================
# CREATE INPUT FILE
# ===========================================
print('rank = ',rank)
if rank == 0:
  print('Writing input file...')

  particleBlockString=""
  particleTypeString=""
  particleRegionString=""
  targetRegionsString=""
  blockIndex = 0
  numMats = len(matsOrig)
  for i in range(numMats):
    numTypes = len(particleTypesPerMat[i])
    types = list(particleTypesPerMat[i])
    regionBlocksStr = ""
    for j in range(numTypes):
      regionBlocksStr += "pb"+str(blockIndex)
      
      if types[j] == 0: # Single point
        particleTypeString+="SinglePoint"
      if types[j] == 1: # Single point with B-splines
        particleTypeString+="SinglePointBSpline"
      if types[j] == 2: # CPDI
        particleTypeString+="CPDI"
      if types[j] == 3: #CPTI
        particleTypeString+="CPTI"
      if types[j] == 4: #CPDI2
        particleTypeString+="CPDI2"
      if types[j] < 0 or types[j] > 4:
        print("Unknown particle type!")
        sys.exit(0)
      
      if j < numTypes-1:
        regionBlocksStr+=", "
      if i < numMats-1:
        particleTypeString+=", "

      blockIndex += 1
    
    particleBlockString += regionBlocksStr
    targetRegionsString += "particles/ParticleRegion"+str(i+1)
    if i < numMats-1:
      particleBlockString+=", "
      targetRegionsString+=", "
    particleRegionString+="""
    <ParticleRegion
        name="ParticleRegion"""+str(i+1)+""""
        meshBody="particles"
        particleBlocks="{ """+ regionBlocksStr + """ }"
        materialList="{ """ + matsOrig[i] + """ }"/>"""

  geosxInputFileName = 'mpm_'+inputFile.replace('pfw_input_',"")+'.xml'
  geosxInputFile = open(geosxInputFileName, 'w')

  geosxInputFileString = """<?xml version="1.0" ?>
<!-- 
srun -p pdebug -n """+str(mNodes)+""" """+geosPath+""" -i """+geosxInputFileName+"""
srun -n """+str(mCores)+""" """+geosPath+""" -i """+geosxInputFileName+"""
-->
<Problem>

  <Mesh>
    <InternalMesh
      name="backgroundGrid"
      elementTypes="{ C3D8 }"
      xCoords=""" + '"{' + str(XMIN) + "," + str(XMAX) + '}"' + """
      yCoords=""" + '"{' + str(YMIN) + "," + str(YMAX) + '}"' + """
      zCoords=""" + '"{' + str(ZMIN) + "," + str(ZMAX) + '}"' + """
      nx="""+'"{'+str(NI)+'}"'+"""
      ny="""+'"{'+str(NJ)+'}"'+"""
      nz="""+'"{'+str(NK)+'}"'+"""
      periodic=""" + '"{ ' + str(int(periodic[0])) + ', ' + str(int(periodic[1])) + ', ' + str(int(periodic[2])) + '}"' + """
      cellBlockNames="{ cb1 }"/>
      
    <ParticleMesh
      name="particles"
      particleFile="""+'"'+particleFileName+'"'+"""
      particleBlockNames="{ """+particleBlockString+""" }"
      particleTypes="{ """+particleTypeString+""" }"/>
  </Mesh>

  <ElementRegions>
    <CellElementRegion
      name="CellRegion1"
      meshBody="backgroundGrid"
      cellBlocks="{ cb1 }"
      materialList="{ null }"/>
  </ElementRegions>

  <ParticleRegions>"""+particleRegionString+"""
  </ParticleRegions>

  <Solvers
    gravityVector="{ 0.0, 0.0, 0.0 }">
    <SolidMechanics_MPM
      name="mpmsolve"
      discretization="FE1"
      targetRegions="{ backgroundGrid/CellRegion1, """+targetRegionsString+""" }"
"""+mpmSolverParameterString+""">
""" + mpmEventsString.replace('\n', '\n      ') + """
    </SolidMechanics_MPM>
  </Solvers>

  <Constitutive>
    <ElasticIsotropic
      name="null"
      defaultDensity="1000"
      defaultBulkModulus="1.0e9"
      defaultShearModulus="1.0e9"/>
    """ + str(materialPropertyString).replace('\n','\n    ') + """
  </Constitutive>

  <Events
    maxTime="""+'"'+str(endTime)+'"'+""">
    <PeriodicEvent
      name="solverApplications"
      target="/Solvers/mpmsolve"/>
    <PeriodicEvent
      name="outputs"
      timeFrequency="""+'"'+str( plotInterval )+'"'+"""
      target="/Outputs/vtkOutput"/>
    <PeriodicEvent
      name="restart"
      timeFrequency=""" + '"' + str( restartInterval ) + '"' + """
      target="/Outputs/restartOutput"/> 
    <HaltEvent
      maxRuntime=""" + '"' + str( maxRestartTime ) + '"' + """/>
  </Events>

  <NumericalMethods>
    <FiniteElements>
      <FiniteElementSpace
        name="FE1"
        order="1"/>
    </FiniteElements>
  </NumericalMethods>

  <Outputs>
    <VTK
      name="vtkOutput"
      format="ascii"/>
    <Restart
      name="restartOutput"/>
  </Outputs>

</Problem>

"""

  geosxInputFile.write(geosxInputFileString)
  geosxInputFile.close()




  # ===========================================
  #  Make a SLURM script to run GEOSX.
  # ===========================================

  # This will run the job.
  # #SBATCH --export=NONE is needed because we are using mpi4py, which modifies
  # the environment, and that somehow gets passed through slurm and messes
  # up the srun command, causing it to hang unless we include this line.

  if lassen:
    slurmScript = """#!/bin/bash
### LSF syntax
#BSUB -nnodes """+str(mNodes)+""" #number of nodes
#BSUB -W """+mWallTimeMinutes+""" #walltime in minutes
#BSUB -G """+mBank+""" #account
#BSUB -J """+str(inputFile+"_RUN")+""" #name of job
#BSUB -q """ + mPartition +  """ #queue to use

echo "Launching jsrun command..."
jsrun -n """+str(mCores)+""" """+geosPath+""" -i """+geosxInputFileName+""" -x """ + str(xpar) + """ -y """ + str(ypar) + """ -z """ + str(zpar) + """
echo "srun command has completed, good bye."
"""
    fileName = timeStamp+"_runGEOSX.sh"
    file = open(fileName, 'w')
    file.write(slurmScript)
    file.close()
  else:  
    slurmScript = """#!/bin/bash
#SBATCH -t """+mWallTime+"""
#SBATCH -N """+str(mNodes)+"""
#SBATCH -A """+mBank+"""
#SBATCH --export=NONE
#SBATCH -p """+ mPartition + """

echo "Launching srun command..."
srun -n """+str(mCores)+""" """+geosPath+""" -i """+geosxInputFileName+""" -x """ + str(xpar) + """ -y """ + str(ypar) + """ -z """ + str(zpar) + """
echo "srun command has completed, good bye."
"""
    fileName = timeStamp+"_runGEOSX.sh"
    file = open(fileName, 'w')
    file.write(slurmScript)
    file.close()

  if mSubmitJobs:

    print('submitting job: '+fileName+' using Popen')
    if lassen:
      output = subprocess.Popen(["bsub", fileName], stdout=subprocess.PIPE).communicate()[0]
    else:
      output = subprocess.Popen(["sbatch", fileName], stdout=subprocess.PIPE).communicate()[0]

    output = str(output, 'UTF8')
    output = output.strip('Submitted batch job ')
    jobID = output.strip()
    print('Submitted job with ID = ',output.strip())

  if ( mSubmitJobs and autoRestart ):
    if lassen:
      print('XXXX  autoRestart not supported with lassen')
    print('Auto restart enabled')
    print('run_check output = ',output.strip())
    slurmScript = """#!/bin/bash
#SBATCH -t 00:02:00
#SBATCH -N 1
#SBATCH -p """+ mPartition +"""
#SBATCH -A """+mBank+"""
#SBATCH --dependency=afterany:"""+jobID+"""
echo "launching pfw_check script..."
python3 pfw_check.py """+inputFile+""" """+jobID+"""
echo "pfw_check script has completed, good bye."

"""
    fileName = timeStamp+"_runCheck.sh"
    file = open(fileName, 'w')
    file.write(slurmScript)
    file.close()    
    call(["sbatch",fileName], cwd=PWD)
