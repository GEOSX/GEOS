

=========================== ================== ======== ================================================================================================================================================== 
Name                        Type               Default  Description                                                                                                                                        
=========================== ================== ======== ================================================================================================================================================== 
capPressureEpsilon          real64             1e-06    Saturation at which the extremum capillary pressure is attained; used to avoid infinite capillary pressure values for saturations close to 0 and 1 
name                        groupName          required A name is required for any non-unique nodes                                                                                                        
phaseCapPressureExponentInv real64_array       {0.5}    Inverse of capillary power law exponent for each phase                                                                                             
phaseCapPressureMultiplier  real64_array       {1}      Entry pressure value for each phase                                                                                                                
phaseMinVolumeFraction      real64_array       {0}      Minimum volume fraction value for each phase                                                                                                       
phaseNames                  groupNameRef_array required List of fluid phases                                                                                                                               
=========================== ================== ======== ================================================================================================================================================== 


