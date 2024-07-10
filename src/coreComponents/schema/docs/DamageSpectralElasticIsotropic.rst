

============================= ============ ======== ==================================================================== 
Name                          Type         Default  Description                                                          
============================= ============ ======== ==================================================================== 
biotCoefficient               real64_array {0}      Biot coefficient                                                     
criticalStrainEnergy          real64       required Critical stress in a 1d tension test                                 
damagePressure                real64       0        A prescribed uniform pressure in the phase-field crack               
defaultBulkModulus            real64       -1       Default Bulk Modulus Parameter                                       
defaultCompressiveStrength    real64       0        Default compressive strength from the uniaxial compression test      
defaultCriticalFractureEnergy real64       required Default critical fracture energy                                     
defaultDeltaCoefficient       real64       -1       Default coefficient in the calculation of the external driving force 
defaultDensity                real64       required Default Material Density                                             
defaultDrainedLinearTEC       real64       0        Default Linear Thermal Expansion Coefficient of the Solid Rock Frame 
defaultPoissonRatio           real64       -1       Default Poisson's Ratio                                              
defaultShearModulus           real64       -1       Default Shear Modulus Parameter                                      
defaultTensileStrength        real64       0        Default tensile strength from the uniaxial tension test              
defaultYoungModulus           real64       -1       Default Young's Modulus                                              
degradationLowerLimit         real64       0        The lower limit of the degradation function                          
extDrivingForceFlag           integer      0        Whether to have external driving force. Can be 0 or 1                
lengthScale                   real64       required Length scale l in the phase-field equation                           
name                          groupName    required A name is required for any non-unique nodes                          
============================= ============ ======== ==================================================================== 


