

======================= ============ ======== =============================================================================================================================================================== 
Name                    Type         Default  Description                                                                                                                                                     
======================= ============ ======== =============================================================================================================================================================== 
computeCFLNumbers       integer      0        Flag to decide whether CFL numbers are computed or not                                                                                                          
computeRegionStatistics integer      1        Flag to decide whether region statistics are computed or not                                                                                                    
flowSolverName          groupNameRef required Name of the flow solver                                                                                                                                         
logLevel                integer      0        Log level                                                                                                                                                       
name                    groupName    required A name is required for any non-unique nodes                                                                                                                     
relpermThreshold        real64       1e-06    Flag to decide whether a phase is considered mobile (when the relperm is above the threshold) or immobile (when the relperm is below the threshold) in metric 2 
writeCSV                integer      0        Write statistics into a CSV file                                                                                                                                
======================= ============ ======== =============================================================================================================================================================== 


