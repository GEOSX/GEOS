.. _AddingNewSolver:

Adding a new Physics Solver
###########################

In this tutorial, you will learn how to construct a new GEOS Physics Solver class.
We will use *LaplaceFEM* solver, computing the solution of the Laplace problem in
a specified material, as a starting point.

.. math::
   :nowrap:

   \begin{eqnarray*}
   D^* \Delta X = f \quad \mbox{in\;} \Omega \\
   X = X^g \quad \mbox{on\;} \Gamma_g
   \end{eqnarray*}

It is advised to read :ref:`XML_and_classes` preliminary to this tutorial.
The goal of this document is to explain how to develop a new solver that solves
Laplace's equation with a constant diffusion coefficient that is specified by users in the XML input.

For readability, member functions in the text will be referenced by their names but their
arguments will be omitted.

*LaplaceFEM* overview
=====================
The *LaplaceFEM* solver can be found in ``./src/coreComponents/physicsSolvers/simplePDE/``.
Let us inspect declarations in ``LaplaceFEM.hpp`` and implementations in ``LaplaceFEM.cpp``
before diving into specifying a new solver class that meets our needs.

Declaration file (reference)
----------------------------
The included header is ``physicsSolvers/simplePDE/LaplaceBaseH1.hpp`` which declares the base class ``LaplaceBaseH1``, shared by all Laplace solvers. Moreover, ``physicsSolver/simplePDE/LaplaceBaseH1.hpp`` includes the following headers:

 - ``common/EnumStrings.hpp`` which includes facilities for enum-string conversion (useful for reading enum values from input);
 - ``physicsSolver/SolverBase.hpp`` which declares the abstraction class shared by all physics solvers.
 - ``managers/FieldSpecification/FieldSpecificationManager.hpp`` which declares a manager used to access and to set field on the discretized domain.

Let us jump forward to the class enum and variable as they contain the data used
specifically in the implementation of *LaplaceFEM*.

class enums and variables (reference)
`````````````````````````````````````
The class exhibits two member variables:

 - ``m_fieldName`` which stores the name of the diffused variable (*e.g.* the temperature) as a `string`;
 - ``m_timeIntegrationOption`` an `enum` value allowing to dispatch with respect to the transient treatment.

``TimeIntegrationOption`` is an `enum` specifying the transient treatment which can be chosen
respectively between *SteadyState* and *ImplicitTransient* depending on whether we are interested in
the transient state.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_TIMEINTOPT
   :end-before: //END_SPHINX_INCLUDE_TIMEINTOPT

In order to register an enumeration type with the Data Repository and have its value read from input,
we must define stream insertion/extraction operators. This is a common task, so GEOS provides
a facility for automating it. Upon including ``common/EnumStrings.hpp``, we can call the following macro
at the namespace scope (in this case, right after the ``LaplaceBaseH1`` class definition is complete):

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_REGENUM
   :end-before: //END_SPHINX_INCLUDE_REGENUM

Once explained the main variables and enum, let us start reading through the different member functions:

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_BEGINCLASS
   :end-before: //END_SPHINX_INCLUDE_BEGINCLASS

Start looking at the class *LaplaceFEM* constructor and destructor declarations
shows the usual `string` ``name`` and `Group*` pointer to ``parent`` that are required
to build the global file-system like structure of GEOS (see :ref:`GroupPar` for details).
It can also be noted that the nullary constructor is deleted on purpose to avoid compiler
automatic generation and user misuse.

The next method ``catalogName()`` is static and returns the key to be added to the *Catalog* for this type of solver
(see :ref:`ObjectCatalogPar` for details). It has to be paired with the following macro in the implementation file.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_REGISTER
   :end-before: //END_SPHINX_INCLUDE_REGISTER

Finally, the member function ``registerDataOnMesh()`` is declared in the ``LaplaceBaseH1`` class as

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_REGISTERDATAONMESH
   :end-before: //END_SPHINX_INCLUDE_REGISTERDATAONMESH

It is used to assign fields onto the discretized mesh object and
will be further discussed in the :ref:`Implementation` section.

The next block consists in solver interface functions. These member functions set up
and specialize every time step from the system matrix assembly to the solver stage.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_SOLVERINTERFACE
   :end-before: //END_SPHINX_INCLUDE_SOLVERINTERFACE

Furthermore, the following functions are inherited from the base class.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_SOLVERINTERFACE
   :end-before: //END_SPHINX_INCLUDE_SOLVERINTERFACE

Eventually, ``applyDirichletBCImplicit()`` is the working specialized member functions called
when ``applyBoundaryConditions()`` is called in this particular class override.

Browsing the base class ``SolverBase``, it can be noted that most of the solver interface functions are called during
either ``SolverBase::linearImplicitStep()`` or ``SolverBase::nonlinearImplicitStep()`` depending on the solver strategy chosen.

Switching to protected members, ``postInputInitialization()`` is a central member function and
will be called by ``Group`` object after input is read from XML entry file.
It will set and dispatch solver variables from the base class ``SolverBase`` to the most derived class.
For ``LaplaceFEM``, it will allow us to set the right time integration scheme based on the XML value
as will be further explored in the next :ref:`Implementation` section.

Let us focus on a ``struct`` that plays an important role: the *viewKeyStruct* structure.

*viewKeyStruct* structure (reference)
`````````````````````````````````````

This embedded instantiated structure is a common pattern shared by all solvers.
It stores ``dataRepository::ViewKey`` type objects that are used as binding data
between the input XML file and the source code.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_VIEWKEY
   :end-before: //END_SPHINX_INCLUDE_VIEWKEY

We can check that in the *LaplaceFEM* companion integratedTest

.. literalinclude:: ../../../../../inputFiles/simplePDE/Laplace_base.xml
   :language: xml
   :start-after: <Solvers>
   :end-before: </Solvers>

In the following section, we will see where this binding takes place.

.. _Implementation:

Implementation File (reference)
-------------------------------
Switching to implementation, we will focus on few implementations, leaving details
to other tutorials. The ``LaplaceFEM`` constructor is implemented as follows.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_CONSTRUCTOR
   :end-before: //END_SPHINX_INCLUDE_CONSTRUCTOR

As we see, it calls the ``LaplaceBaseH1`` constructor, that is implemented as follows.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_CONSTRUCTOR
   :end-before: //END_SPHINX_INCLUDE_CONSTRUCTOR

Checking out the constructor, we can see that the use of a ``registerWrapper<T>(...)``
allows us to register the key value from the `enum` ``viewKeyStruct`` defining them as:

 - ``InputFlags::OPTIONAL`` if they are optional and can be provided;
 - ``InputFlags::REQUIRED`` if they are required and will throw error if not;

and their associated descriptions for auto-generated docs.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceBaseH1.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_REGISTERDATAONMESH
   :end-before: //END_SPHINX_INCLUDE_REGISTERDATAONMESH

``registerDataOnMesh()`` is browsing all subgroups in the mesh ``Group`` object and
for all nodes in the sub group:

 - register the observed field under the chosen ``m_fieldName`` key;
 - apply a default value;
 - set the output verbosity level (here ``PlotLevel::LEVEL_0``);
 - set the field associated description for auto generated docs.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/simplePDE/LaplaceFEM.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_ASSEMBLY
   :end-before: //END_SPHINX_INCLUDE_ASSEMBLY

``assembleSystem()`` will be our core focus as we want to change the diffusion coefficient from its
hard coded value to a XML read user-defined value. One can see that this method is in charge of constructing
in a parallel fashion the FEM system matrix. Bringing ``nodeManager`` and ``ElementRegionManager`` from domain local
``MeshLevel`` object together with ``FiniteElementDiscretizationManager`` from the ``NumericalMethodManager``, it uses
nodes embedded loops on degrees of freedom in a local index embedded loops to fill a matrix and a rhs container.

As we spotted the place to change in a code to get a user-defined diffusion coefficient into the game, let us jump
to writing our new *LaplaceDiffFEM* solver.

.. note::

  We might want to remove final keyword from ``postInputInitialization()`` as it will prevent you from overriding it.

Start doing your own Physic solver
==================================
As we will extend *LaplaceFEM* capabilities, we will derive publicly from it.

Declaration File
----------------

As there is only few places where we have to change, the whole declaration file is reported below and
commented afterwards.

.. code-block:: c++

  #include "physicsSolvers/simplePDE/LaplaceFEM.hpp"

  namespace geos
  {

  class LaplaceDiffFEM : public LaplaceFEM
  {
  public:

    LaplaceDiffFEM() = delete;

    LaplaceDiffFEM( const string& name,
                    Group * const parent );

    virtual ~LaplaceDiffFEM() override;

    static string catalogName() { return "LaplaceDiffFEM"; }

    virtual void
    assembleSystem( real64 const time,
                    real64 const dt,
                    DomainPartition * const domain,
                    DofManager const & dofManager,
                    ParallelMatrix & matrix,
                    ParallelVector & rhs ) override;


    struct viewKeyStruct : public LaplaceFEM::viewKeyStruct
    {
      dataRepository::ViewKey diffusionCoeff = { "diffusionCoeff" };
    } laplaceDiffFEMViewKeys;

    protected:
    virtual void postInputInitialization() override final;

  private:
    real64 m_diffusion;

  };


We intend to have a user-defined diffusion coefficient, we then need a `real64` class variable ``m_diffusion``
to store it.

Consistently with *LaplaceFEM*, we will also delete the nullary constructor and declare a constructor with the same arguments for
forwarding to `Group` master class. Another mandatory step is to override the static ``CatalogName()`` method to properly
register any data from the new solver class.

Then as mentioned in :ref:`Implementation`, the diffusion coefficient is used when assembling the matrix coefficient. Hence
we will have to override the ``assembleSystem()`` function as detailed below.

Moreover, if we want to introduce a new binding between the input XML and the code we will have to work on the three
``struct viewKeyStruct`` , ``postInputInitialization()`` and the constructor.

Our new solver ``viewKeyStruct`` will have its own structure inheriting from the *LaplaceFEM* one to have the ``timeIntegrationOption``
and ``fieldName`` field. It will also create a ``diffusionCoeff`` field to be bound to the user defined homogeneous coefficient on one hand
and to our ``m_diffusion`` class variable on the other.


Implementation File
-------------------
As we have seen in :ref:`Implementation`, the first place where to implement a new register from XML input is
in the constructor. The ``diffusionCoeff`` entry we have defined in the ``laplaceDiffFEMViewKeys``
will then be asked as a required input. If not provided, the error thrown will ask for it described asked
an "input uniform diffusion coefficient for the Laplace equation".

.. code-block:: c++

  LaplaceDiffFEM::LaplaceDiffFEM( const string& name,
                                  Group * const parent ):
  LaplaceFEM( name, parent ), m_diffusion(0.0)
  {
    registerWrapper<string>(laplaceDiffFEMViewKeys.diffusionCoeff.Key()).
      setInputFlag(InputFlags::REQUIRED).
      setDescription("input uniform diffusion coeff for the laplace equation");
  }

Another important spot for binding the value of the XML read parameter to our ``m_diffusion`` is in ``postInputInitialization()``.

.. code-block:: c++

  void LaplaceDiffFEM::postInputInitialization()
  {
    LaplaceFEM::postInputInitialization();

    string sDiffCoeff = this->getReference<string>(laplaceDiffFEMViewKeys.diffusionCoeff);
    this->m_diffusion = std::stof(sDiffCoeff);
  }

Now that we have required, read and bind the user-defined diffusion value to a variable, we can use it in the construction of our
matrix into the overridden ``assembleSystem()``.

.. code-block:: c++
  :emphasize-lines: 16-18

  // begin element loop, skipping ghost elements
  for( localIndex k=0 ; k<elementSubRegion->size() ; ++k )
  {
    if(elemGhostRank[k] < 0)
    {
      element_rhs = 0.0;
      element_matrix = 0.0;
      for( localIndex q=0 ; q<n_q_points ; ++q)
      {
        for( localIndex a=0 ; a<numNodesPerElement ; ++a)
        {
          elemDofIndex[a] = dofIndex[ elemNodes( k, a ) ];

          for( localIndex b=0 ; b<numNodesPerElement ; ++b)
          {
            element_matrix(a,b) += detJ[k][q] *
                                   m_diffusion *
                                 + Dot( dNdX[k][q][a], dNdX[k][q][b] );
          }

        }
      }
      matrix.add( elemDofIndex, elemDofIndex, element_matrix );
      rhs.add( elemDofIndex, element_rhs );
    }
  }

This completes the implementation of our new solver *LaplaceDiffFEM*.

Nonetheless, the compiler should complain that ``m_fieldName`` is privately as inherited from *LaplaceFEM*. One should then either promote ``m_fieldName`` to protected
or add a getter in *LaplaceFEM* class to correct the error. The getter option has been chosen and the fix in our solver is then:

.. code-block:: c++

  array1d<globalIndex> const & dofIndex =
    nodeManager->getReference< array1d<globalIndex> >( dofManager.getKey( getFieldName() ) );


Note: For consistency do not forget to change LaplaceFEM to LaplaceDiffFEM in the guards comments

Last steps
==========

After assembling both declarations and implementations for our new solver, the final steps go as:

 - add declarations to parent CMakeLists.txt (here add to ``physicsSolvers_headers`` );
 - add implementations to parent CMakeLists.txt (here add to ``physicsSolvers_sources``);
 - check that Doxygen comments are properly set in our solver class;
 - uncrustify it to match the code style by going to the build folder and running the command: make uncrustify_style;
 - write unit tests for each new features in the solver class;
 - write an integratedTests for the solver class.
