.. _LvArray:

LvArray
##################################################

Use in GEOS
=============================

``LvArray`` containers are used in GEOS as primary storage mechanism for mesh topology, field data and any other type of "large" data sets (i.e. ones that scale with the size of the problem).
When allocating a new field, using one of ``LvArray`` containers is mandatory if the data is meant to be used in any computational kernel.
The file ``common/DataTypes.hpp`` provides shorthand aliases for commonly used containers:

.. literalinclude:: ../../../../coreComponents/common/DataTypes.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: //END_SPHINX_INCLUDE_00


LvArray documentation
=====================

Please refer to the full `LvArray documentation <https://lvarray.readthedocs.io/en/latest/>`_ for details on each of the classes.
