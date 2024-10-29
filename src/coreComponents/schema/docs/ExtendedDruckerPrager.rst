

============================ ========= ======== ==================================================================== 
Name                         Type      Default  Description                                                          
============================ ========= ======== ==================================================================== 
defaultBulkModulus           real64    -1       Default Bulk Modulus Parameter                                       
defaultCohesion              real64    0        Initial cohesion                                                     
defaultDensity               real64    required Default Material Density                                             
defaultDilationRatio         real64    1        Dilation ratio [0,1] (ratio = tan dilationAngle / tan frictionAngle) 
defaultDrainedLinearTEC      real64    0        Default Linear Thermal Expansion Coefficient of the Solid Rock Frame 
defaultHardening             real64    0        Hardening parameter (hardening rate is faster for smaller values)    
defaultInitialFrictionAngle  real64    30       Initial friction angle (degrees)                                     
defaultPoissonRatio          real64    -1       Default Poisson's Ratio                                              
defaultResidualFrictionAngle real64    30       Residual friction angle (degrees)                                    
defaultShearModulus          real64    -1       Default Shear Modulus Parameter                                      
defaultYoungModulus          real64    -1       Default Young's Modulus                                              
name                         groupName required A name is required for any non-unique nodes                          
============================ ========= ======== ==================================================================== 


