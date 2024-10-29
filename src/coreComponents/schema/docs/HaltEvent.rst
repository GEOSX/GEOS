

==================== ============ ======== ===================================================================================================================== 
Name                 Type         Default  Description                                                                                                           
==================== ============ ======== ===================================================================================================================== 
beginTime            real64       0        Start time of this event.                                                                                             
endTime              real64       1e+100   End time of this event.                                                                                               
finalDtStretch       real64       0.001    Allow the final dt request for this event to grow by this percentage to match the endTime exactly.                    
forceDt              real64       -1       While active, this event will request this timestep value (ignoring any children/targets requests).                   
logLevel             integer      0        Log level                                                                                                             
maxEventDt           real64       -1       While active, this event will request a timestep <= this value (depending upon any child/target requests).            
maxRuntime           real64       required The maximum allowable runtime for the job.                                                                            
name                 groupName    required A name is required for any non-unique nodes                                                                           
target               groupNameRef          Name of the object to be executed when the event criteria are met.                                                    
targetExactStartStop integer      1        If this option is set, the event will reduce its timestep requests to match any specified beginTime/endTimes exactly. 
HaltEvent            node                  :ref:`XML_HaltEvent`                                                                                                  
PeriodicEvent        node                  :ref:`XML_PeriodicEvent`                                                                                              
SoloEvent            node                  :ref:`XML_SoloEvent`                                                                                                  
==================== ============ ======== ===================================================================================================================== 


