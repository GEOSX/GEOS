

==================== ============ ======== ================================================================================================================================================================================= 
Name                 Type         Default  Description                                                                                                                                                                       
==================== ============ ======== ================================================================================================================================================================================= 
beginTime            real64       0        Start time of this event.                                                                                                                                                         
cycleFrequency       integer      1        Event application frequency (cycle, default)                                                                                                                                      
endTime              real64       1e+100   End time of this event.                                                                                                                                                           
finalDtStretch       real64       0.001    Allow the final dt request for this event to grow by this percentage to match the endTime exactly.                                                                                
forceDt              real64       -1       While active, this event will request this timestep value (ignoring any children/targets requests).                                                                               
function             groupNameRef          Name of an optional function to evaluate when the time/cycle criteria are met.If the result is greater than the specified eventThreshold, the function will continue to execute.  
logLevel             integer      0        Log level                                                                                                                                                                         
maxEventDt           real64       -1       While active, this event will request a timestep <= this value (depending upon any child/target requests).                                                                        
name                 groupName    required A name is required for any non-unique nodes                                                                                                                                       
object               groupNameRef          If the optional function requires an object as an input, specify its path here.                                                                                                   
set                  groupNameRef          If the optional function is applied to an object, specify the setname to evaluate (default = everything).                                                                         
stat                 integer      0        If the optional function is applied to an object, specify the statistic to compare to the eventThreshold.The current options include: min, avg, and max.                          
target               groupNameRef          Name of the object to be executed when the event criteria are met.                                                                                                                
targetExactStartStop integer      1        If this option is set, the event will reduce its timestep requests to match any specified beginTime/endTimes exactly.                                                             
targetExactTimestep  integer      1        If this option is set, the event will reduce its timestep requests to match the specified timeFrequency perfectly: dt_request = min(dt_request, t_last + time_frequency - time)). 
threshold            real64       0        If the optional function is used, the event will execute if the value returned by the function exceeds this threshold.                                                            
timeFrequency        real64       -1       Event application frequency (time).  Note: if this value is specified, it will override any cycle-based behavior.                                                                 
HaltEvent            node                  :ref:`XML_HaltEvent`                                                                                                                                                              
PeriodicEvent        node                  :ref:`XML_PeriodicEvent`                                                                                                                                                          
SoloEvent            node                  :ref:`XML_SoloEvent`                                                                                                                                                              
==================== ============ ======== ================================================================================================================================================================================= 


