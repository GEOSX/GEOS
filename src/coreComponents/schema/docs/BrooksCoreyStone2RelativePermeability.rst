

======================= ================== ======== ========================================================================================================================================================== 
Name                    Type               Default  Description                                                                                                                                                
======================= ================== ======== ========================================================================================================================================================== 
gasOilRelPermExponent   real64_array       {1}      | Rel perm power law exponent for the pair (gas phase, oil phase) at residual water saturation                                                               
                                                    | The expected format is "{ gasExp, oilExp }", in that order                                                                                                 
gasOilRelPermMaxValue   real64_array       {0}      | Maximum rel perm value for the pair (gas phase, oil phase) at residual water saturation                                                                    
                                                    | The expected format is "{ gasMax, oilMax }", in that order                                                                                                 
name                    groupName          required A name is required for any non-unique nodes                                                                                                                
phaseMinVolumeFraction  real64_array       {0}      Minimum volume fraction value for each phase                                                                                                               
phaseNames              groupNameRef_array required List of fluid phases                                                                                                                                       
waterOilRelPermExponent real64_array       {1}      | Rel perm power law exponent for the pair (water phase, oil phase) at residual gas saturation                                                               
                                                    | The expected format is "{ waterExp, oilExp }", in that order                                                                                               
waterOilRelPermMaxValue real64_array       {0}      | Maximum rel perm value for the pair (water phase, oil phase) at residual gas saturation                                                                    
                                                    | The expected format is "{ waterMax, oilMax }", in that order                                                                                               
======================= ================== ======== ========================================================================================================================================================== 


