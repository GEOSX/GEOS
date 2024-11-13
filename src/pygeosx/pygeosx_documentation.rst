.. _pygeosxInterface:

:mod:`pygeosx` --- GEOS in Python
==================================
GEOS can be manipulated and executed through a Python script.

High-level control of GEOS is managed through the top-level ``pygeosx`` functions,
like ``initialize`` and ``run``. GEOS's data can be manipulated by getting
`pylvarray <https://lvarray.readthedocs.io/en/latest/python/index.html>`_ views of LvArray objects living in GEOS's data repository.
These ``pylvarray`` views are fetched by calling ``Wrapper.value()`` after getting
a ``Wrapper`` out of the data repository.

.. warning:: 
  The ``pygeosx`` module provides plenty of opportunities to crash Python. 
  See the Segmentation Faults section below.

Only Python 3 is supported.

Module Functions
----------------

.. py:function:: pygeosx.initialize(rank, args)

  Initialize GEOS for the first time, with a rank and command-line arguments.

  This function should only be called **once**. To reinitialize, use the ``reinit`` function.

  Generally the rank is obtained from the `mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_ 
  module and the arguments are obtained from ``sys.argv``.

  Returns a ``Group`` representing the ``ProblemManager`` instance.

.. py:function:: pygeosx.reinit(args)

  Reinitialize GEOS with a new set of command-line arguments.

  Returns a ``Group`` representing the ``ProblemManager`` instance.

.. py:function:: pygeosx.apply_initial_conditions()

  Apply the initial conditions.

.. py:function:: pygeosx.finalize()

  Finalize GEOS. After this no calls into pygeosx or to MPI are allowed.

.. py:function:: pygeosx.run()

  Enter the GEOS event loop.

  Runs until hitting a breakpoint defined in the input deck, or until the simulation
  is complete.

  Returns one of the state constants defined below.

GEOS State
-----------

.. py:data:: pygeosx.UNINITIALIZED

.. py:data:: pygeosx.INITIALIZED

.. py:data:: pygeosx.READY_TO_RUN

  This state indicates that GEOS still has time steps left to run.

.. py:data:: pygeosx.COMPLETED

  This state indicates that GEOS has completed the current simulation.

Module Classes
--------------

.. py:class:: pygeosx.Group

  Python interface to geos::dataRepository::Group.

  Used to get access to other groups, and ultimately to get wrappers
  and convert them into Python views of C++ objects.

  .. py:method:: groups()

    Return a list of the subgroups.

  .. py:method:: wrappers()

    Return a list of the wrappers.

  .. py:method:: get_group(path)
    :noindex:

  .. py:method:: get_group(path, default)

    Return the ``Group`` at the relative path ``path``; ``default`` is optional.
    If no group exists and ``default`` is not given, raise a ``ValueError``;
    otherwise return ``default``.

  .. py:method:: get_wrapper(path)
    :noindex:

  .. py:method:: get_wrapper(path, default)

    Return the ``Wrapper`` at the relative path ``path``; ``default`` is optional.
    If no Wrapper exists and ``default`` is not given, raise a ``ValueError``;
    otherwise return ``default``.

  .. py:method:: register(callback)

    Register a callback on the physics solver.

    The callback should take two arguments: the CRSMatrix and the array.

    Raise ``TypeError`` if the group is not the Physics solver.

.. py:class:: pygeosx.Wrapper

  Python interface to geos::dataRepository::WrapperBase.

  Wraps a generic C++ object. Use ``repr`` to get a description of the type.

  .. py:method:: value()

    Return a view of the wrapped value, or ``None`` if it cannot be exported to Python.

    A breakdown of the possible return types:

    - Instance of a pylvarray class
      If the wrapped type is one of the LvArray types that have a Python wrapper type.

    - 1D numpy.ndarray
      If the wrapped type is a numeric constant. The returned array is a shallow copy and
      has a single entry.

    - str
      If the wrapped type is a std::string this returns a copy of the string.

    - list of str
      If the wrapped type is a `LvArray::Array< std::string, 1, ... >` or a `std::vector< std::string >`.
      This is a copy.

    - None
      If the wrapped type is not covered by any of the above.


Segmentation Faults
-------------------
Improper use of this module and associated programs can easily cause Python to crash. 
There are two main causes of crashes. Both can be avoided by following some
general guidelines.

Stale Numpy Views
^^^^^^^^^^^^^^^^^
The ``pylvarray`` classes (which may be returned from ``Wrapper.value()``) 
provide various ways to get Numpy views of
their data. However, those views are only valid as long as the 
LvArray object's buffer is not reallocated. The buffer may be reallocated
by invoking methods (the ones that require 
the ``pylarray.RESIZEABLE`` permission) or by calls into ``pygeosx``. 
It is strongly recommended that you do not keep Numpy views of 
LvArray objects around after calls to ``pygeosx``.

.. code:: python

  my_array = pygeosx.get_wrapper("path").value()
  view = my_array.to_numpy()
  my_array.resize(1000)
  print(view) # segfault

Destroyed LvArray C++ objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As mentioned earlier, the classes defined in this
module cannot be created in Python; ``pygeosx`` must create an LvArray 
object in C++, then create a
``pylvarray`` view of it. However, the Python view will only be
valid as long as the underlying LvArray C++ object is kept around. If 
that is destroyed, the Python object will be left holding an invalid
pointer and subsequent attempts to use the Python object will cause 
undefined behavior. Unfortunately, ``pygeosx`` may destroy LvArray
objects without warning. It is therefore strongly recommended that you do
not keep ``pylvarray`` objects around after calls to ``pygeosx``. The following
code snippet, for instance, could segfault:

.. code:: python

  my_array = pygeosx.get_wrapper("path").value()
  pygeosx.run()
  view = my_array.to_numpy() # segfault


.. toctree::
   :maxdepth: 2
   :caption: Contents:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

