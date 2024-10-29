

=============== ================== =========== =============================================================================== 
Name            Type               Default     Description                                                                     
=============== ================== =========== =============================================================================== 
childDirectory  string                         Child directory path                                                            
filename        string             TimeHistory The filename to which to write time history output.                             
format          string             hdf         The output file format for time history output.                                 
logLevel        integer            0           Log level                                                                       
name            groupName          required    A name is required for any non-unique nodes                                     
parallelThreads integer            1           Number of plot files.                                                           
sources         groupNameRef_array required    A list of collectors from which to collect and output time history information. 
=============== ================== =========== =============================================================================== 


