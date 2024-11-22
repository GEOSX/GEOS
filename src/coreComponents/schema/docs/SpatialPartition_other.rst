

================== ==================================================================================================================================================== ================================================================ 
Name               Type                                                                                                                                                 Description                                                      
================== ==================================================================================================================================================== ================================================================ 
blockSize          real64_array                                                                                                                                         Length of partition dimensions (excluding ghost objects).        
contactGhostMax    real64_array                                                                                                                                         Ghost position max.                                              
contactGhostMin    real64_array                                                                                                                                         Ghost position min.                                              
gridMax            real64_array                                                                                                                                         Maximum extent of problem dimensions (excluding ghost objects).  
gridMin            real64_array                                                                                                                                         Minimum extent of problem dimensions (excluding ghost objects).  
gridSize           real64_array                                                                                                                                         Total length of problem dimensions (excluding ghost objects).    
max                real64_array                                                                                                                                         Maximum extent of partition dimensions (excluding ghost objects) 
min                real64_array                                                                                                                                         Minimum extent of partition dimensions (excluding ghost objects) 
partitionLocations LvArray_Array< LvArray_Array< double, 1, camp_int_seq< long, 0l >, int, LvArray_ChaiBuffer >, 1, camp_int_seq< long, 0l >, int, LvArray_ChaiBuffer > Locations of partition boundaries                                
================== ==================================================================================================================================================== ================================================================ 


