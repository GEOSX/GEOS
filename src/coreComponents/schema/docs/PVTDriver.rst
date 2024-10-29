

====================== ============ ======== ===================================================================== 
Name                   Type         Default  Description                                                           
====================== ============ ======== ===================================================================== 
baseline               path         none     Baseline file                                                         
feedComposition        real64_array required Feed composition array [mol fraction]                                 
fluid                  groupNameRef required Fluid to test                                                         
logLevel               integer      0        Log level                                                             
name                   groupName    required A name is required for any non-unique nodes                           
output                 string       none     Output file                                                           
outputCompressibility  integer      0        Flag to indicate that the total compressibility should be output      
outputMassDensity      integer      0        Flag to indicate that the mass density of each phase should be output 
outputPhaseComposition integer      0        Flag to indicate that phase compositions should be output             
pressureControl        groupNameRef required Function controlling pressure time history                            
steps                  integer      required Number of load steps to take                                          
temperatureControl     groupNameRef required Function controlling temperature time history                         
====================== ============ ======== ===================================================================== 


