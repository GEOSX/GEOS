

=============== ==================================== ======== =========================================================================== 
Name            Type                                 Default  Description                                                                 
=============== ==================================== ======== =========================================================================== 
coordinateFiles path_array                           {}       List of coordinate file names for ND Table                                  
coordinates     real64_array                         {0}      Coordinates inputs for 1D tables                                            
inputVarNames   groupNameRef_array                   {}       Name of fields are input to function.                                       
interpolation   geos_TableFunction_InterpolationType linear   | Interpolation method. Valid options:                                        
                                                              | * linear                                                                    
                                                              | * nearest                                                                   
                                                              | * upper                                                                     
                                                              | * lower                                                                     
name            groupName                            required A name is required for any non-unique nodes                                 
values          real64_array                         {0}      Values for 1D tables                                                        
voxelFile       path                                          Voxel file name for ND Table                                                
=============== ==================================== ======== =========================================================================== 


