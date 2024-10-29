

========================= ================================================================= ====================================================================== 
Name                      Type                                                              Description                                                            
========================= ================================================================= ====================================================================== 
domainBoundaryIndicator   integer_array                                                     (no description available)                                             
ghostRank                 integer_array                                                     (no description available)                                             
globalToLocalMap          geos_mapBase<long long, int, std_integral_constant<bool, false> > (no description available)                                             
isExternal                integer_array                                                     (no description available)                                             
localToGlobalMap          globalIndex_array                                                 Array that contains a map from localIndex to globalIndex.              
location                  real64_array2d                                                    For each perforation, physical location (x,y,z coordinates)            
numPerforationsGlobal     globalIndex                                                       (no description available)                                             
reservoirElementIndex     integer_array                                                     For each perforation, element index of the perforated element          
reservoirElementRegion    integer_array                                                     For each perforation, elementRegion index of the perforated element    
reservoirElementSubregion integer_array                                                     For each perforation, elementSubRegion index of the perforated element 
wellElementIndex          integer_array                                                     For each perforation, index of the well element                        
wellSkinFactor            real64_array                                                      For each perforation, well skin factor                                 
wellTransmissibility      real64_array                                                      For each perforation, well transmissibility                            
neighborData              node                                                              :ref:`DATASTRUCTURE_neighborData`                                      
sets                      node                                                              :ref:`DATASTRUCTURE_sets`                                              
========================= ================================================================= ====================================================================== 


