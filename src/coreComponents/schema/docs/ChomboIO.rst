

================== ========= =================== =============================================================================================== 
Name               Type      Default             Description                                                                                     
================== ========= =================== =============================================================================================== 
beginCycle         real64    required            Cycle at which the coupling will commence.                                                      
childDirectory     string                        Child directory path                                                                            
inputPath          string    /INVALID_INPUT_PATH Path at which the chombo to geosx file will be written.                                         
name               groupName required            A name is required for any non-unique nodes                                                     
outputPath         string    required            Path at which the geosx to chombo file will be written.                                         
parallelThreads    integer   1                   Number of plot files.                                                                           
useChomboPressures integer   0                   True iff geosx should use the pressures chombo writes out.                                      
waitForInput       integer   required            True iff geosx should wait for chombo to write out a file. When true the inputPath must be set. 
================== ========= =================== =============================================================================================== 


