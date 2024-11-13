.. _WrapperDoc:


Wrapper
=======

This class encapsulates an object for storage in a ``Group`` and provides an interface for performing some common operations on that object.

Description
-----------

In the filesystem analogy, a Wrapper may be thought of as a file that stores actual data.
Each ``Wrapper`` belong to a single ``Group`` much like a file belongs to a filesystem directory.
In general, more than one wrapper in the tree may refer to the same wrapped object, just like symlinks in the file system may refer to the same file.
However, only one wrapper should be *owning* the data (see below).

In the XML input file, ``Wrapper`` correspond to attribute of an XML element representing the containing ``Group``.
See :ref:`XML_and_classes` for the relationship between XML input files and Data Repository.

``Wrapper<T>`` is templated on the type of object it encapsulates, thus providing strong type safety when retrieving the objects.
As each ``Wrapper`` class instantiation will be a distinct type, Wrapper derives from a non-templated ``WrapperBase`` class that defines a common interface.
``WrapperBase`` is the type of pointer that is stored in the ``MappedVector`` container within a ``Group``.

``WrapperBase`` provides several interface functions that delegate the work to the wrapped object if it supports the corresponding method signature.
This allows a collection of heterogeneous wrappers (i.e. over different types) to be treated uniformly.
Examples include:

* ``size()``
* ``resize(newSize)``
* ``reserve(newCapacity)``
* ``capacity()``
* ``move(LvArray::MemorySpace)``

A ``Wrapper`` may be *owning* or *non-owning*, depending on how it's constructed.
An *owning* ``Wrapper`` will typically either take a previously allocated object via ``std::unique_ptr<T>`` or no pointer at all and itself allocate the object.
It will delete the wrapped object when destroyed.
A *non-owning* ``Wrapper`` may be used to register with the data repository objects that are not directly heap-allocated, for example data members of other objects.
It will take a raw pointer as input and not delete the wrapped object when destroyed.

Attributes
----------

Each instance of ``Wrapper`` has a set of attributes that control its function in the data repository.
These attributes are:

* **InputFlags**

  A strongly typed enum that defines the relationship between the Wrapper and the XML input.
  Possible values are:

  +--------------+-----------------------------------------------------------------------------+
  | Value        | Explanation                                                                 |
  +==============+=============================================================================+
  | ``FALSE``    | Data is not read from XML input (default).                                  |
  +--------------+-----------------------------------------------------------------------------+
  | ``OPTIONAL`` | Data is read from XML if an attribute matching Wrapper's name is found.     |
  +--------------+-----------------------------------------------------------------------------+
  | ``REQUIRED`` | Data is read from XML and an error is raised if the attribute is not found. |
  +--------------+-----------------------------------------------------------------------------+

  Other values of ``InputFlags`` enumeration are reserved for ``Group`` objects.

.. note::
   A runtime error will occur when attempting to read from XML a wrapped type ``T`` that does not have ``operator>>`` defined.

* **RestartFlags**

  Enumeration that describes how the Wrapper interacts with restart files.

  +--------------------+---------------------------------------------------------------+
  | Value              | Explanation                                                   |
  +====================+===============================================================+
  | ``NO_WRITE``       | Data is not written into restart files.                       |
  +--------------------+---------------------------------------------------------------+
  | ``WRITE``          | Data is written into restart files but not read upon restart. |
  +--------------------+---------------------------------------------------------------+
  | ``WRITE_AND_READ`` | Data is both written and read upon restart (default).         |
  +--------------------+---------------------------------------------------------------+

.. note::
   A runtime error will occur when attempting to write a wrapped type ``T`` that does not support buffer packing.
   Therefore, when registering custom types (i.e. not a basic C++ type or an :ref:`LvArray` container) we recommend setting the flag to ``NO_WRITE``.
   A future documentation topic will explain how to extend buffer packing capabilities to custom user-defined types.

* **PlotLevel**

  Enumeration that describes how the Wrapper interacts with plot (visualization) files.

  +-------------+-------------------------------------------------------------------+
  | Value       | Explanation                                                       |
  +=============+===================================================================+
  | ``LEVEL_0`` | Data always written to plot files.                                |
  +-------------+-------------------------------------------------------------------+
  | ``LEVEL_1`` | Data written to plot when ``plotLevel``>=1 is specified in input. |
  +-------------+-------------------------------------------------------------------+
  | ``LEVEL_2`` | Data written to plot when ``plotLevel``>=2 is specified in input. |
  +-------------+-------------------------------------------------------------------+
  | ``LEVEL_3`` | Data written to plot when ``plotLevel``>=3 is specified in input. |
  +-------------+-------------------------------------------------------------------+
  | ``NOPLOT``  | Data never written to plot files.                                 |
  +-------------+-------------------------------------------------------------------+

.. note::
   Only data stored in :ref:`LvArray`'s ``Array<T>`` containers is currently written into plot files.

Default Values
--------------

``Wrapper`` supports setting a default value for its wrapped object.
The default value is used if a wrapper with ``InputFlags::OPTIONAL`` attribute does not match an attribute in the input file.
For :ref:`LvArray` containers it is also used as a default value for new elements upon resizing the container.

Default value can be set via one of the following two methods:

* ``setDefaultValue`` sets the default value but does not affect the actual value stored in the wrapper.
* ``setApplyDefaultValue`` sets the default value *and* applies it to the stored value.

.. note::
   A runtime error is raised if a default value is not set for a wrapper with ``InputFlags::OPTIONAL`` attribute.

The type ``DefaultValue<T>`` is used to store the default value for the wrapper.

.. todo::
   ``DefaultValue`` is actually not a type but an alias for another internal struct.
   As such, it cannot currently be specialized for a user's custom type.

Doxygen API documentation
-------------------------

`Wrapper API <../../../doxygen_output/html/classgeos_1_1data_repository_1_1_wrapper.html>`_
