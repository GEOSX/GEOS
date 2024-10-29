

=============================== =================================== ========== ===================================================================== 
Name                            Type                                Default    Description                                                           
=============================== =================================== ========== ===================================================================== 
maxPermeability                 real64                              1          Max. permeability can be reached.                                     
name                            groupName                           required   A name is required for any non-unique nodes                           
pressureDependenceConstants     R1Tensor                            required   Pressure dependence coefficients for each permeability component.     
pressureModelType               geos_constitutive_PressureModelType Hyperbolic Type of the pressure dependence model.                                
referencePermeabilityComponents R1Tensor                            required   Reference xx, yy and zz components of a diagonal permeability tensor. 
referencePressure               real64                              required   Reference pressure for the pressure permeability model                
=============================== =================================== ========== ===================================================================== 


