

============== ================== ======== ===================================================================================================================================================================================== 
Name           Type               Default  Description                                                                                                                                                                           
============== ================== ======== ===================================================================================================================================================================================== 
flowSolverName groupNameRef       required Name of the flow solver                                                                                                                                                               
fluxNames      groupNameRef_array {*}      Name(s) array of the SourceFlux(s) for which we want the statistics. Use "*" to target all SourceFlux.                                                                                
logLevel       integer            0        | Log level                                                                                                                                                                             
                                           | - Log Level 1 outputs the sum of all SourceFlux(s) produced rate & mass,                                                                                                              
                                           | - Log Level 2 details values for each SourceFlux,                                                                                                                                     
                                           | - Log Level 3 details values for each region.                                                                                                                                         
name           groupName          required A name is required for any non-unique nodes                                                                                                                                           
writeCSV       integer            0        Write statistics into a CSV file                                                                                                                                                      
============== ================== ======== ===================================================================================================================================================================================== 


