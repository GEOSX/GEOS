.. _Group:


Group
=====

``dataRepository::Group`` serves as a base class for most objects in GEOS.
In GEOS, the ``Group`` may be thought of as an analogy to the file folder in a hierachical filesystem-like structure.
As such, a ``Group`` is used as a container class that holds a collection of other Groups, or sub-Groups,
a pointer to the parent of the Group, and a collection of Wrappers.
The ``Group`` also defines a general capability to create and traverse/access the objects in the hierarchy.
The Wrappers contained in a Group may be of arbitrary type, but in the case of an LvArray object, a Group
size and capacity may be translated down to array, thereby keeping a collection of wrapped Array objects
that will resize in unison with the Group.
Each group has a string "name" that defines its key in the parent Group.
This key string must be unique in the scope of the parent Group.

Implementation Details
----------------------
Some noteworthy implementation details inside the declaration of ``dataRepository::Group`` are:

.. literalinclude:: ../Group.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: //END_SPHINX_INCLUDE_00
   
* In the GEOS repository, the ``keyType`` is specified to be a ``string`` for all  collection objects, 
  while the ``indexType`` is specified to be a ``localIndex``.
  The types are set in the ``common/DataTypes.hpp`` file, but are typically a ``string`` and a 
  ``std::ptrdiff_t`` respectively.
  
.. literalinclude:: ../Group.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_01
   :end-before: //END_SPHINX_INCLUDE_01
   
* The ``subGroupMap`` and ``wrapperMap`` aliases represent the type of container that the collection of 
  sub-``Group`` s and ``Wrapper`` s are stored in for each ``Group``.
  These container types are template specializations of the ``MappedVector`` class, which store a pointer to 
  a type, and provides functionality for a key or index based lookup. 
  More details may be found in the documentation for ``MappedVector``.
  
.. literalinclude:: ../Group.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_02
   :end-before: //END_SPHINX_INCLUDE_02
   
* The ``m_parent`` member is a pointer to the ``Group`` that contains the current ``Group`` as part of its
  collection of sub-``Group`` s.

  .. warning::
    The existence of the non-const ``m_parent`` gives the current ``Group`` access to alter the parent Group.
    Special care should be taken to avoid using this access whenever possible. 
    Remember...with great power comes great responsibility.
* The ``m_wrappers`` member is the collection of Wrappers contained in the current ``Group``.
* The ``m_subGroups`` member is the collection of ``Group`` s contained in the current ``Group``.
* The ``m_size`` and ``m_capacity`` members are used to set the size and capacity of any objects contained
  in the ``m_wrappers`` collection that have been specified to be set by their owning ``Group``.
  This is typically only useful for Array types and is implemented within the ``WrapperBase`` object.
* The ``m_name`` member is the key of this Group in the collection of ``m_parent->m_subGroups``.
  This key is unique in the scope of ``m_parent``, so some is required when constructing the hierarchy.
  
Interface Functions
-------------------
The public interface for ``dataRepository::Group`` provides functionality for constructing a hierarchy, 
and traversing that hierarchy, as well as accessing the contents of objects stored in the ``Wrapper`` 
containers stored within a ``Group``.

Adding New Groups
^^^^^^^^^^^^^^^^^
To add new sub-``Group`` s there are several ``registerGroup`` functions that add a new ``Group`` under
the calling ``Group`` scope.
A listing of these functions is provided:

.. literalinclude:: ../Group.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_REGISTER_GROUP
   :end-before: //END_SPHINX_INCLUDE_REGISTER_GROUP

These functions all take in a ``name`` for the new ``Group``, which will be used as the key when trying to 
access the ``Group`` in the future.
Some variants create a new ``Group``, while some variants take in an existing ``Group`` . 
The template argument is to specify the actaul type of the ``Group`` as it it is most likely a type that 
derives from ``Group`` that is we would like to create in the repository.
Please see the doxygen documentation for a detailed description of each option.

Getting Groups
^^^^^^^^^^^^^^
The collection of functions to retrieve a ``Group`` and their descriptions are taken from source and shown 
here:

.. literalinclude:: ../Group.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_GET_GROUP
   :end-before: //END_SPHINX_INCLUDE_GET_GROUP


Register Wrappers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../Group.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_REGISTER_WRAPPER
   :end-before: //END_SPHINX_INCLUDE_REGISTER_WRAPPER


Getting Wrappers/Wrapped Objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../Group.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_GET_WRAPPER
   :end-before: //END_SPHINX_INCLUDE_GET_WRAPPER
   
Looping Interface
^^^^^^^^^^^^^^^^^

.. literalinclude:: ../Group.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_LOOP_INTERFACE
   :end-before: //END_SPHINX_INCLUDE_LOOP_INTERFACE


Doxygen API documentation
-------------------------

`Group API <../../../doxygen_output/html/classgeos_1_1data_repository_1_1_group.html>`_
