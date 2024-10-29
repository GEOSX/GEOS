

=========================== ================== ======== ============================================================================================================================================== 
Name                        Type               Default  Description                                                                                                                                    
=========================== ================== ======== ============================================================================================================================================== 
capPressureEpsilon          real64             1e-06    Wetting-phase saturation at which the max cap. pressure is attained; used to avoid infinite cap. pressure values for saturations close to zero 
name                        groupName          required A name is required for any non-unique nodes                                                                                                    
phaseCapPressureExponentInv real64_array       {2}      Inverse of capillary power law exponent for each phase                                                                                         
phaseEntryPressure          real64_array       {1}      Entry pressure value for each phase                                                                                                            
phaseMinVolumeFraction      real64_array       {0}      Minimum volume fraction value for each phase                                                                                                   
phaseNames                  groupNameRef_array required List of fluid phases                                                                                                                           
=========================== ================== ======== ============================================================================================================================================== 


