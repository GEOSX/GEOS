

============= ======================== ======== ===================================================================== 
Name          Type                     Default  Description                                                           
============= ======================== ======== ===================================================================== 
axialControl  groupNameRef             required Function controlling axial stress or strain (depending on test mode)  
baseline      path                     none     Baseline file                                                         
initialStress real64                   required Initial stress (scalar used to set an isotropic stress state)         
logLevel      integer                  0        Log level                                                             
material      groupNameRef             required Solid material to test                                                
mode          geos_TriaxialDriver_Mode required Test mode [stressControl, strainControl, mixedControl]                
name          groupName                required A name is required for any non-unique nodes                           
output        string                   none     Output file                                                           
radialControl groupNameRef             required Function controlling radial stress or strain (depending on test mode) 
steps         integer                  required Number of load steps to take                                          
============= ======================== ======== ===================================================================== 


