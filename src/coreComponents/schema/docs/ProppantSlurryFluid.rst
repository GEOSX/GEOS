

========================= ============ ======== ================================================ 
Name                      Type         Default  Description                                      
========================= ============ ======== ================================================ 
componentNames            string_array {}       List of fluid component names                    
compressibility           real64       0        Fluid compressibility                            
defaultComponentDensity   real64_array {0}      Default value for the component density.         
defaultComponentViscosity real64_array {0}      Default value for the component viscosity.       
defaultCompressibility    real64_array {0}      Default value for the component compressibility. 
flowBehaviorIndex         real64_array {0}      Flow behavior index                              
flowConsistencyIndex      real64_array {0}      Flow consistency index                           
maxProppantConcentration  real64       0.6      Maximum proppant concentration                   
name                      groupName    required A name is required for any non-unique nodes      
referenceDensity          real64       1000     Reference fluid density                          
referencePressure         real64       100000   Reference pressure                               
referenceProppantDensity  real64       1400     Reference proppant density                       
referenceViscosity        real64       0.001    Reference fluid viscosity                        
========================= ============ ======== ================================================ 


