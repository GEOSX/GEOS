

=========================== ================== ======== ======================================================================================================================================================================================== 
Name                        Type               Default  Description                                                                                                                                                                              
=========================== ================== ======== ======================================================================================================================================================================================== 
childDirectory              string                      Child directory path                                                                                                                                                                     
fieldNames                  groupNameRef_array {}       Names of the fields to output. If this attribute is specified, GEOSX outputs all (and only) the fields specified by the user, regardless of their plotLevel                              
name                        groupName          required A name is required for any non-unique nodes                                                                                                                                              
onlyPlotSpecifiedFieldNames integer            0        If this flag is equal to 1, then we only plot the fields listed in `fieldNames`. Otherwise, we plot all the fields with the required `plotLevel`, plus the fields listed in `fieldNames` 
parallelThreads             integer            1        Number of plot files.                                                                                                                                                                    
plotFileRoot                string             plot     (no description available)                                                                                                                                                               
plotLevel                   integer            1        (no description available)                                                                                                                                                               
writeCellElementMesh        integer            1        (no description available)                                                                                                                                                               
writeEdgeMesh               integer            0        (no description available)                                                                                                                                                               
writeFEMFaces               integer            0        (no description available)                                                                                                                                                               
writeFaceElementMesh        integer            1        (no description available)                                                                                                                                                               
=========================== ================== ======== ======================================================================================================================================================================================== 


