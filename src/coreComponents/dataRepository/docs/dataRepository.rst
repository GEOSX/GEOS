.. _dataRepository:

###############
Data Repository
###############

The GEOS "Data Repository" is intended to provide the building blocks for the code structure within GEOS.
The "Data Repository" provides a general capability to store arbitrary data and objects in a hierarchical
structure, similar to a standard file system. The "Wrapper" object is a generic container for any object,
and provides a standard interface for accessing an object and performing standard operations. 
The "Group" object is a container for "Wrapper" and other "Group" objects.

The components/classes of the data structure that a developer should have some knowledge of are:

.. toctree::
   :maxdepth: 2

   /coreComponents/dataRepository/docs/Group
   /coreComponents/dataRepository/docs/Wrapper
   /coreComponents/dataRepository/docs/ObjectCatalog
   /coreComponents/dataRepository/docs/MappedVector
