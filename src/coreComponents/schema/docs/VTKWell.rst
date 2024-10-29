

===================== ========= ======== ==================================================================================================== 
Name                  Type      Default  Description                                                                                          
===================== ========= ======== ==================================================================================================== 
file                  path      required Path to the well file                                                                                
minElementLength      real64    0.001    Minimum length of a well element, computed as (segment length / number of elements per segment ) [m] 
minSegmentLength      real64    0.01     Minimum length of a well segment [m]                                                                 
name                  groupName required A name is required for any non-unique nodes                                                          
numElementsPerSegment integer   required Number of well elements per polyline segment                                                         
radius                real64    required Radius of the well [m]                                                                               
wellControlsName      string    required Name of the set of constraints associated with this well                                             
wellRegionName        string    required Name of the well element region                                                                      
Perforation           node               :ref:`XML_Perforation`                                                                               
===================== ========= ======== ==================================================================================================== 


