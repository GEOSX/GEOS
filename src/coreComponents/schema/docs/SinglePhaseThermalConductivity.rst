

===================================== ========= ======== =============================================================================================================== 
Name                                  Type      Default  Description                                                                                                     
===================================== ========= ======== =============================================================================================================== 
defaultThermalConductivityComponents  R1Tensor  required xx, yy, and zz diagonal components of the default thermal conductivity tensor [J/(s.m.K)]                       
name                                  groupName required A name is required for any non-unique nodes                                                                     
referenceTemperature                  real64    0        The reference temperature at which the conductivity components are equal to the default values                  
thermalConductivityGradientComponents R1Tensor  {0,0,0}  xx, yy, and zz diagonal components of the thermal conductivity gradient tensor w.r.t. temperature [J/(s.m.K^2)] 
===================================== ========= ======== =============================================================================================================== 


