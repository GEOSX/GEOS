

============================================ ================== ======== ==================================================================================================== 
Name                                         Type               Default  Description                                                                                          
============================================ ================== ======== ==================================================================================================== 
bcApplicationTableName                       groupNameRef                Name of table that specifies the on/off application of the boundary condition.                       
beginTime                                    real64             -1e+99   Time at which the boundary condition will start being applied.                                       
componentFractionVsElevationTableNames       groupNameRef_array {}       Names of the tables specifying the (component fraction vs elevation) relationship for each component 
componentNames                               string_array       {}       Names of the fluid components                                                                        
datumElevation                               real64             required Datum elevation [m]                                                                                  
datumPressure                                real64             required Datum pressure [Pa]                                                                                  
direction                                    R1Tensor           {0,0,0}  Direction to apply boundary condition to.                                                            
elevationIncrementInHydrostaticPressureTable real64             0.6096   Elevation increment [m] in the hydrostatic pressure table constructed internally                     
endTime                                      real64             1e+99    Time at which the boundary condition will stop being applied.                                        
equilibrationTolerance                       real64             0.001    Tolerance in the fixed-point iteration scheme used for hydrostatic initialization                    
functionName                                 groupNameRef                Name of function that specifies variation of the boundary condition.                                 
initialPhaseName                             groupNameRef                Name of the phase initially saturating the reservoir                                                 
logLevel                                     integer            0        Log level                                                                                            
maxNumberOfEquilibrationIterations           integer            5        Maximum number of equilibration iterations                                                           
name                                         groupName          required A name is required for any non-unique nodes                                                          
objectPath                                   groupNameRef                Path to the target field                                                                             
scale                                        real64             0        Scale factor for value of the boundary condition.                                                    
temperatureVsElevationTableName              groupNameRef                Name of the table specifying the (temperature [K] vs elevation) relationship                         
============================================ ================== ======== ==================================================================================================== 


