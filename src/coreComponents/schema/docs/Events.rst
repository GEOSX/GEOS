

================ ================================== ============ ======================================================================== 
Name             Type                               Default      Description                                                              
================ ================================== ============ ======================================================================== 
logLevel         integer                            0            Log level                                                                
maxCycle         integer                            2147483647   Maximum simulation cycle for the global event loop. Disabled by default. 
maxTime          real64                             1.79769e+308 Maximum simulation time for the global event loop. Disabled by default.  
minTime          real64                             0            Start simulation time for the global event loop.                         
timeOutputFormat geos_EventManager_TimeOutputFormat seconds      Format of the time in the GEOS log.                                      
HaltEvent        node                                            :ref:`XML_HaltEvent`                                                     
PeriodicEvent    node                                            :ref:`XML_PeriodicEvent`                                                 
SoloEvent        node                                            :ref:`XML_SoloEvent`                                                     
================ ================================== ============ ======================================================================== 


