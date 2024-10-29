

================================= ================== ======== ================================================================================== 
Name                              Type               Default  Description                                                                        
================================= ================== ======== ================================================================================== 
name                              groupName          required A name is required for any non-unique nodes                                        
phaseNames                        groupNameRef_array required List of fluid phases                                                               
phaseThermalConductivity          real64_array       required Phase thermal conductivity [W/(m.K)]                                               
rockThermalConductivityComponents R1Tensor           required xx, yy, and zz components of a diagonal rock thermal conductivity tensor [W/(m.K)] 
================================= ================== ======== ================================================================================== 


