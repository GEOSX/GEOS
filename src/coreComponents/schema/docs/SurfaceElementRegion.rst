

=============== ============================================== ================= ================================================================================= 
Name            Type                                           Default           Description                                                                       
=============== ============================================== ================= ================================================================================= 
defaultAperture real64                                         required          The default aperture of newly formed surface elements.                            
faceBlock       groupNameRef                                   FractureSubRegion The name of the face block in the mesh, or the embedded surface.                  
materialList    groupNameRef_array                             required          List of materials present in this region                                          
meshBody        groupNameRef                                                     Mesh body that contains this region                                               
name            groupName                                      required          A name is required for any non-unique nodes                                       
subRegionType   geos_SurfaceElementRegion_SurfaceSubRegionType faceElement       Type of surface element subregion. Valid options: {faceElement, embeddedElement}. 
=============== ============================================== ================= ================================================================================= 


