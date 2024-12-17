.. _MappedVectorDoc:


Mapped Vector
=============

A superposition of a contiguous and an associative container.

Description
-----------

The container stores pointers to objects (which are themselves heap-allocated).
Each element may be optionally *owned* by the container, in which case it will be deleted upon removal or container destruction.
The pointers are stored in a contiguous memory allocation, and thus are accessible through an integral index lookup.
In addition, there is a map that provides a key lookup capability to the container if that is the preferred interface.

The container template has four type parameters:

* ``T`` is the object type pointed to by container entries

* ``T_PTR`` is a pointer-to-``T`` type which must be either ``T *`` (default) or ``std::unique_ptr<T>``

* ``KEY_TYPE`` is the type of key used in associative lookup

* ``INDEX_TYPE`` is the type used in index lookup

Element access
--------------

``MappedVector`` provides three main types of data access using ``[]`` operator:

* **Index lookup** is the fastest way of element random access if the ordinal index is known.

* **Key lookup** is similar to key lookup of any associative container and incurs similar cost.

* **KeyIndex lookup** uses a special type, ``KeyIndex``, that contains both a key and an index.
  Initially the index is unknown and the key is used for the lookup.
  The ``KeyIndex`` is modified during lookup, storing the index located.
  If the user persists the ``KeyIndex`` object, they may reuse it in subsequent accesses and get the benefit of direct index access.

In addition to these, an STL-conformant iterator interface is available via ``begin()`` and ``end()`` methods.
The type iterated over is a key-pointer pair (provided as `value_type` alias).

Doxygen API documentation
-------------------------

`MappedVector  API <../../../doxygen_output/html/classgeos_1_1_mapped_vector.html>`_
