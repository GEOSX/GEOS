

=========================== ======================================= ======== ======================================================================================== 
Name                        Type                                    Default  Description                                                                              
=========================== ======================================= ======== ======================================================================================== 
collisionAlpha              real64                                  1.27     Collision alpha coefficient                                                              
collisionBeta               real64                                  1.5      Collision beta coefficient                                                               
fluidViscosity              real64                                  0.001    Fluid viscosity                                                                          
hinderedSettlingCoefficient real64                                  5.9      Hindered settling coefficient                                                            
isCollisionalSlip           integer                                 0        Whether the collisional component of the slip velocity is considered                     
maxProppantConcentration    real64                                  0.6      Max proppant concentration                                                               
name                        groupName                               required A name is required for any non-unique nodes                                              
particleSettlingModel       geos_constitutive_ParticleSettlingModel required | Particle settling velocity model. Valid options:                                         
                                                                             | * Stokes                                                                                 
                                                                             | * Intermediate                                                                           
                                                                             | * Turbulence                                                                             
proppantDensity             real64                                  1400     Proppant density                                                                         
proppantDiameter            real64                                  0.0002   Proppant diameter                                                                        
slipConcentration           real64                                  0.1      Slip concentration                                                                       
sphericity                  real64                                  1        Sphericity                                                                               
=========================== ======================================= ======== ======================================================================================== 


