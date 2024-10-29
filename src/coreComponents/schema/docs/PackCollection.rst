

====================== ================== ======== =========================================================================================== 
Name                   Type               Default  Description                                                                                 
====================== ================== ======== =========================================================================================== 
disableCoordCollection integer            0        Whether or not to create coordinate meta-collectors if collected objects are mesh objects.  
fieldName              groupNameRef       required The name of the (packable) field associated with the specified object to retrieve data from 
name                   groupName          required A name is required for any non-unique nodes                                                 
objectPath             groupNameRef       required The name of the object from which to retrieve field values.                                 
onlyOnSetChange        integer            0        Whether or not to only collect when the collected sets of indices change in any way.        
setNames               groupNameRef_array {}       The set(s) for which to retrieve data.                                                      
====================== ================== ======== =========================================================================================== 


