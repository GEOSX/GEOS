

======================== ============================= ======== ================================================================================== 
Name                     Type                          Default  Description                                                                        
======================== ============================= ======== ================================================================================== 
adiosConfig              string                                 Path to the adios configuration file when using the catalyst-adios implementation. 
childDirectory           string                                 Child directory path                                                               
fullFieldChannelName     string                                 Name to give to the channel passing the full field data.                           
implementation           string                                 Name of the catalyst implementation to use.                                        
implementationPath       string                                 Path to the catalyst the implementation to use.                                    
name                     string                        required A name is required for any non-unique nodes                                        
outputFullQuadratureData integer                       0        If true writes out data associated with every quadrature point.                    
parallelThreads          integer                       1        Number of plot files.                                                              
plotLevel                geos_dataRepository_PlotLevel 1        Determines which fields to write.                                                  
scripts                  string                        required Column separated paths to the catalyst scripts.                                    
======================== ============================= ======== ================================================================================== 


