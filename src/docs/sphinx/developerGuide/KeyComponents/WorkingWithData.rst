.. _WorkingWithData:

#####################################
Working with data in GEOS
#####################################

In ``GEOS``, data is typically registered in the :ref:`dataRepository`.
This allows for the writing/reading of data to/from restart and plot files.
Any object that derives from :ref:`Group` may have data registered on it
through the methods described in :ref:`Group`.
Similarly, accessing data from outside the scope of an object is possible
through one of the various ``Group::get()``.
Of course, for temporary data that does not need to persist between cycles,
or across physics packages, you may simply define member or local variable which
will not be registered with the :ref:`dataRepository`.

Working with data on the Mesh objects
=====================================
The mesh objects in ``GEOS`` such as the ``FaceManager`` or ``NodeManager``,
are derived from ``ObjectManagerBase``, which in turn derives from :ref:`Group`.
The important distinction is that ``ObjectManagerBase`` contains various members
that are useful when defining mesh object managers.
When considering data that is attached to a mesh object, we group the data into
two categories:

* `Intrinsic <https://www.merriam-webster.com/dictionary/intrinsic>`_ data is
  data that is required to describe the object.
  For instance, to define a ``Node``, the ``NodeManager`` contains an array of
  positions corresponding to each ``Node`` it contains.
  Thus the ``ReferencePosition`` is ``Intrinsic`` data.
  ``Intrinsic`` data is almost always a member of the mesh object, and is
  registered on the mesh object in the constructor of mesh object itself.
* Field data (or `Extrinsic <https://www.merriam-webster.com/dictionary/extrinsic>`_ data) is
  data that is not required to define the object.
  For instance, a physics package
  may request that a ``Velocity`` value be stored on the nodes.
  Appropriately the data will be registered on the ``NodeManager``.
  However, this data is not required to define a ``Node``, and is viewed as
  ``Fields`` or ``Extrinsic``.
  ``Field`` data is never a member of the mesh object, and is typically
  registered on the mesh object outside of the definition of the mesh object
  (i.e. from a physics solver).


Registering Intrinsic data on a Mesh Object
-------------------------------------------
As mentioned above, ``Intrinsic`` data is typically a member of the mesh object,
and is registered in the constructor of the mesh Object.
Taking the ``NodeManager`` and the ``referencePosition`` as an example, we
point out that the reference position is actually a member in the
``NodeManager``.

.. literalinclude:: ../../../../coreComponents/mesh/NodeManager.hpp
   :language: c++
   :start-after: //START_SPHINX_REFPOS
   :end-before: //END_SPHINX_REFPOS

This member is registered in the constructor for the ``NodeManager``.

.. literalinclude:: ../../../../coreComponents/mesh/NodeManager.cpp
   :language: c++
   :start-after: //START_SPHINX_REFPOS_REG
   :end-before: //END_SPHINX_REFPOS_REG

Finally in order to access this data, the ``NodeManager`` provides explicit
accessors.

.. literalinclude:: ../../../../coreComponents/mesh/NodeManager.hpp
   :language: c++
   :start-after: //START_SPHINX_REFPOS_ACCESS
   :end-before: //END_SPHINX_REFPOS_ACCESS

Thus the interface for ``Intrinsic`` data is set by the object that it is a part
of, and the developer may only access the data through the accesssors from
outside of the mesh object class scope.

Registering Field data on a Mesh Object
---------------------------------------
To register ``Field`` data, there are many ways a developer may proceed.
We will use the example of registering a ``totalDisplacement`` on the ``NodeManager``
from the ``SolidMechanics`` solver.
The most general approach is to define a string key and call one of the
`Group::registerWrapper() <../../../doxygen_output/html/classgeos_1_1data_repository_1_1_group.html#a741c3b5728fc47b33fbaad6c4f124991>`_
functions from ``PhysicsSolverBase::registerDataOnMesh()``.
Then when you want to use the data, you can call ``Group::getReference()``.
For example this would look something like:

.. code-block:: c++

    void SolidMechanicsLagrangianFEM::registerDataOnMesh( Group * const MeshBodies )
    {
      for( auto & mesh : MeshBodies->GetSubGroups() )
      {
        NodeManager & nodes = mesh.second->groupCast< MeshBody * >()->getMeshLevel( 0 ).getNodeManager();

        nodes.registerWrapper< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( keys::totalDisplacement ).
          setPlotLevel( PlotLevel::LEVEL_0 ).
          setRegisteringObjects( this->getName()).
          setDescription( "An array that holds the total displacements on the nodes." ).
          reference().resizeDimension< 1 >( 3 );
      }
    }

and

.. code-block:: c++

    arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & u = nodes.getReference< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( keys::totalDisplacement );
    ... do something with u

This approach is flexible and extendible, but is potentially error prone due to
its verbosity and lack of information centralization.
Therefore we also provide a more controlled/uniform method by which to register
and extract commonly used data on the mesh.
The ``trait approach`` requires the definition of a ``traits struct`` for each
data object that will be supported.
To apply the ``trait approach`` to the example use case shown above, there
should be the following definition somewhere in a header file:

.. code-block:: c++

    namespace fields
    {
    struct totalDisplacement
    {
      static constexpr auto key = "totalDisplacement";
      using DataType = real64;
      using Type = array2d< DataType, nodes::TOTAL_DISPLACEMENT_PERM >;
      static constexpr DataType defaultValue = 0;
      static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0;

      /// Description of the data associated with this trait.
      static constexpr auto description = "An array that holds the total displacements on the nodes.";
    };
    }

Also note that you should use the ``DECLARE_FIELD`` C++ macro that will perform this tedious task for you.
Then the registration is simplified as follows:

.. code-block:: c++

    void SolidMechanicsLagrangianFEM::registerDataOnMesh( Group * const MeshBodies )
    {
      for( auto & mesh : MeshBodies->GetSubGroups() )
      {
        NodeManager & nodes = mesh.second->groupCast< MeshBody * >()->getMeshLevel( 0 ).getNodeManager();
        nodes.registerField< fields::totalDisplacement >( this->getName() ).resizeDimension< 1 >( 3 );
      }
    }

And to extract the data, the call would be:

.. code-block:: c++

    arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & u = nodes.getField< fields::totalDisplacement >();
    ... do something with u

The end result of the ``trait approach`` to this example is that the developer
has defined a standard specification for ``totalDisplacement``, which may be
used uniformly across the code.
