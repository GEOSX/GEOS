

============= ================== ======== ========================================================================== 
Name          Type               Default  Description                                                                
============= ================== ======== ========================================================================== 
expression    string                      Composite math expression                                                  
functionNames string_array       {}       List of source functions. The order must match the variableNames argument. 
inputVarNames groupNameRef_array {}       Name of fields are input to function.                                      
name          groupName          required A name is required for any non-unique nodes                                
variableNames groupNameRef_array {}       List of variables in expression                                            
============= ================== ======== ========================================================================== 


