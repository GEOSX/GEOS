

================= ================== ======== ======================================================================================================= 
Name              Type               Default  Description                                                                                             
================= ================== ======== ======================================================================================================= 
cellBlockNames    groupNameRef_array required Names of each mesh block                                                                                
elementTypes      string_array       required Element types of each mesh block                                                                        
name              groupName          required A name is required for any non-unique nodes                                                             
nx                integer_array      required Number of elements in the x-direction within each mesh block                                            
ny                integer_array      required Number of elements in the y-direction within each mesh block                                            
nz                integer_array      required Number of elements in the z-direction within each mesh block                                            
positionTolerance real64             1e-10    A position tolerance to verify if a node belong to a nodeset                                            
trianglePattern   integer            0        Pattern by which to decompose the hex mesh into wedges                                                  
xBias             real64_array       {1}      Bias of element sizes in the x-direction within each mesh block (dx_left=(1+b)*L/N, dx_right=(1-b)*L/N) 
xCoords           real64_array       required x-coordinates of each mesh block vertex                                                                 
yBias             real64_array       {1}      Bias of element sizes in the y-direction within each mesh block (dy_left=(1+b)*L/N, dx_right=(1-b)*L/N) 
yCoords           real64_array       required y-coordinates of each mesh block vertex                                                                 
zBias             real64_array       {1}      Bias of element sizes in the z-direction within each mesh block (dz_left=(1+b)*L/N, dz_right=(1-b)*L/N) 
zCoords           real64_array       required z-coordinates of each mesh block vertex                                                                 
InternalWell      node                        :ref:`XML_InternalWell`                                                                                 
VTKWell           node                        :ref:`XML_VTKWell`                                                                                      
================= ================== ======== ======================================================================================================= 


