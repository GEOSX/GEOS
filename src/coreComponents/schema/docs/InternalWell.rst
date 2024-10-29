

===================== =================== ======== ==================================================================================================== 
Name                  Type                Default  Description                                                                                          
===================== =================== ======== ==================================================================================================== 
logLevel              integer             0        Log level                                                                                            
minElementLength      real64              0.001    Minimum length of a well element, computed as (segment length / number of elements per segment ) [m] 
minSegmentLength      real64              0.01     Minimum length of a well segment [m]                                                                 
name                  groupName           required A name is required for any non-unique nodes                                                          
numElementsPerSegment integer             required Number of well elements per polyline segment                                                         
polylineNodeCoords    real64_array2d      required Physical coordinates of the well polyline nodes                                                      
polylineSegmentConn   globalIndex_array2d required Connectivity of the polyline segments                                                                
radius                real64              required Radius of the well [m]                                                                               
wellControlsName      string              required Name of the set of constraints associated with this well                                             
wellRegionName        string              required Name of the well element region                                                                      
Perforation           node                         :ref:`XML_Perforation`                                                                               
===================== =================== ======== ==================================================================================================== 


