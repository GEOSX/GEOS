

======================= =========================================== ======== =================================================================================== 
Name                    Type                                        Default  Description                                                                         
======================= =========================================== ======== =================================================================================== 
compressibility         real64                                      0        Fluid compressibility                                                               
defaultDensity          real64                                      required Default value for density.                                                          
defaultViscosity        real64                                      required Default value for viscosity.                                                        
densityModelType        geos_constitutive_ExponentApproximationType linear   | Type of density model. Valid options:                                               
                                                                             | * exponential                                                                       
                                                                             | * linear                                                                            
                                                                             | * quadratic                                                                         
internalEnergyModelType geos_constitutive_ExponentApproximationType linear   | Type of internal energy model. Valid options:                                       
                                                                             | * exponential                                                                       
                                                                             | * linear                                                                            
                                                                             | * quadratic                                                                         
name                    groupName                                   required A name is required for any non-unique nodes                                         
referenceDensity        real64                                      1000     Reference fluid density                                                             
referenceInternalEnergy real64                                      0.001    Reference fluid internal energy                                                     
referencePressure       real64                                      0        Reference pressure                                                                  
referenceTemperature    real64                                      0        Reference temperature                                                               
referenceViscosity      real64                                      0.001    Reference fluid viscosity                                                           
specificHeatCapacity    real64                                      0        Fluid heat capacity. Unit: J/kg/K                                                   
thermalExpansionCoeff   real64                                      0        Fluid thermal expansion coefficient. Unit: 1/K                                      
viscosibility           real64                                      0        Fluid viscosity exponential coefficient                                             
viscosityModelType      geos_constitutive_ExponentApproximationType linear   | Type of viscosity model. Valid options:                                             
                                                                             | * exponential                                                                       
                                                                             | * linear                                                                            
                                                                             | * quadratic                                                                         
======================= =========================================== ======== =================================================================================== 


