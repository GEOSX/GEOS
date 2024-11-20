###############################################################################
DoF Manager
###############################################################################

This will contains a description of the DoF manager in GEOS.

Brief description
========================

The main aim of the Degrees-of-Freedom (DoF) Manager class is to handle all
degrees of freedom associated with fields that exist on mesh elements, faces, edges and nodes.
It creates a map between local mesh objects and global DoF indices.
Additionally, DofManager simplifies construction of system matrix sparsity patterns.

Key concepts are locations and connectors.
Locations, that can be elements, faces, edges or nodes, represent where the DoF is assigned.
For example, a DoF for pressure in a two-point flux approximation will be on a cell (i.e. element), while a displacement DoF for structural equations will be on a node.
The counterparts of locations are connectors, that are the geometrical entities
that link together different DoFs that create the sparsity pattern.
Connectors can be elements, faces, edges, nodes or none.
Using the same example as before, connectors will be faces and cells, respectively.
The case of a mass matrix, where every element is linked only to itself, is an example when there are no connectors, i.e. these have to be set to none.

DoFs located on a mesh object are owned by the same rank that owns the object in parallel mesh partitioning.
Two types of DoF numbering are supported, with the difference only showing in parallel runs of multi-field problems.

  * Initially, each field is assigned an independent DoF numbering that starts at 0 and is contiguous across all MPI ranks.
    Within each rank, locally owned DoFs are numbered sequentially across mesh locations, and within each mesh location (e.g. node) - sequentially according to component number.
    With this numbering, sparsity patterns can be constructed for individual sub-matrices that represent diagonal/off-diagonal blocks of the global coupled system matrix.

  * After all fields have been declared, the user can call ``DofManager::reorderByRank()``, which constructs a globally contiguous DoF numbering across all fields.
    Specifically, all DoFs owned by rank 0 are numbered field-by-field starting from 0, then those on rank 1, etc.
    This makes global system sparsity pattern compatible with linear algebra packages that only support contiguous matrix rows on each rank.
    At this point, coupled system matrix sparsity pattern can be constructed.

Thus, each instance of ``DofManager`` only supports one type of numbering.
If both types are required, the user is advised to maintain two separate instances of ``DofManager``.


``DofManager`` allocates a separate "DOF index" array for each field on the mesh.
It is an array of global indices, where each value represents the first DoF index for that field and location (or equivalently, the row and column offset of that location's equations and variables for the field in the matrix).
For example, if index array for a field with 3 components contains the value N, global DoF numbers for that location will be N, N+1, N+2.
DoF on ghosted locations have the same indices as on the owning rank.
The array is stored under a generated key, which can be queried from the DoF manager, and is typically used in system assembly.

Methods
========================

The main methods of ``DoF Manager`` are:

* ``setDomain``: sets the domain containing mesh bodies to operate on
  ``domain`` identifies the global domain

 .. code-block:: c

  void setDomain( DomainPartition * const domain );

* ``addField``: creates a new set of DoF, labeled ``field``, with specific
  ``location``.
  Default number of ``components`` is ``1``, like for pressure in flux.
  Default ``regions`` is the empty string, meaning all domain.

 .. code-block:: c

  void addField( string const & fieldName,
                 Location const location,
                 localIndex const components,
                 arrayView1d< string const > const & regions );

* ``addCoupling``: creates a coupling between two fields (``rowField`` and
  ``colField``) according to a given ``connectivity`` in the regions defined by ``regions``.
  Both fields (row and column) must have already been defined on the regions where is required the coupling among them.
  Default value for ``regions`` is the whole intersection between the regions where the first and the second fields are defined.
  This method also creates the coupling between ``colField`` and ``rowField``, i.e. the transpose of the rectangular sparsity pattern.
  This default behaviour can be disabled by passing ``symmetric = false``.

 .. code-block:: c

  void addCoupling( string const & rowField,
                    string const & colField,
                    Connectivity const connectivity,
                    arrayView1d< string const > const & regions,
                    bool const symmetric );

* ``reorderByRank``: finish populating field and coupling information and apply DoF
  re-numbering

 .. code-block:: c

  void reorderByRank();

* ``getKey``: returns the "key" associated with the field, that can be used to access the index array on the mesh object manager corresponding to field's location.

 .. code-block:: c

  string const & getKey( string const & fieldName );

* ``clear``: removes all fields, releases memory and re-opens the DofManager

 .. code-block:: c

  void clear();

* ``setSparsityPattern``: populates the sparsity for the given
  ``rowField`` and ``colField`` into ``matrix``.
  Closes the matrix if ``closePattern`` is ``true``.

 .. code-block:: c

  void setSparsityPattern( MATRIX & matrix,
                           string const & rowField,
                           string const & colField,
                           bool closePattern = true) const;

* ``setSparsityPattern``: populates the sparsity for the full system matrix into ``matrix``.
  Closes the matrix if ``closePattern`` is ``true``.

 .. code-block:: c

  void setSparsityPattern( MATRIX & matrix,
                           bool closePattern = true ) const;

* ``numGlobalDofs``: returns the total number of DoFs across all processors for
  the specified name ``field`` (if given) or all fields (if empty).

 .. code-block:: c

  globalIndex numGlobalDofs( string const & field = "" ) const;

* ``numLocalDofs``: returns the number of DoFs on this process for the
  specified name ``field`` (if given) or all fields (if empty).

 .. code-block:: c

  localIndex numLocalDofs( string const & field = "" ) const;

* ``printFieldInfo``: prints a short summary of declared fields and coupling to the output stream ``os``.

 .. code-block:: c

  void printFieldInfo( std::ostream & os = std::cout ) const;

Example
=======

Here we show how the sparsity pattern is computed for a simple 2D quadrilateral mesh with 6 elements.
Unknowns are pressure, located on the element center, and displacements (*x* and *y* components), located on the nodes.
For fluxes, a two-point flux approximation (TPFA) is used.
The representation of the sparsity pattern of the :math:`\mathsf{C_L}` matrix (connectors/locations) for the simple mesh, shown in :numref:`meshDofManagerFig`, is
reported in :numref:`CLDofManagerFig`.
It can be noticed that the two unknowns for the displacements *x* and *y* are grouped together.
Elements are the connectivity for DoF on nodes (Finite Element Method for displacements) and on elements (pressures).
Faces are the connectivity for DoF on elements (Finite Volume Method for pressure), being the flux computation based on the pressure on the two adjacent elements.

.. _meshDofManagerFig:
.. figure:: /coreComponents/linearAlgebra/docs/images/mesh2D.svg
   :align: center
   :width: 250
   :figclass: align-center

   Small 2D quadrilateral mesh used for this examples.
   Nodes are label with black numbers, elements with light gray numbers and
   faces with italic dark gray numbers.

.. _CLDofManagerFig:
.. figure:: /coreComponents/linearAlgebra/docs/images/CL.svg
   :align: center
   :width: 500
   :figclass: align-center

   Sparsity pattern of the binary matrix connections/locations.

The global sparsity pattern, shown in :numref:`patternDofManagerFig`, is obtained through the symbolic multiplication of the transpose of the matrix :math:`\mathsf{C_L}` and the matrix itself, i.e. :math:`\mathsf{P = C_L^T C_L}`.

.. _patternDofManagerFig:
.. figure:: /coreComponents/linearAlgebra/docs/images/pattern.svg
   :align: center
   :width: 400
   :figclass: align-center

   Sparsity pattern of the global matrix, where red and green entries are related to the displacement field and to the pressure field, respectively.
   Blue entries represent coupling blocks.

Real mesh and patterns
======================

Now we build the pattern of the Jacobian matrix for a simple 3D mesh, shown in
:numref:`meshCubeDofManagerFig`. Fields are:

- displacement (location: node, connectivity: element) defined on the blue, orange and red regions;
- pressure (location: element, connectivity: face) defined on the green, orange and red regions;
- mass matrix (location: element, connectivity: element) defined on the green region only.

Moreover, following coupling are imposed:

- displacement-pressure (connectivity: element) on the orange region only;
- pressure-mass matrix and transpose (connectivity: element) everywhere it is
  possibile.

.. _meshCubeDofManagerFig:
.. figure:: /coreComponents/linearAlgebra/docs/images/meshCube3D.svg
   :align: center
   :width: 400
   :figclass: align-center

   Real mesh used to compute the Jacobian pattern.

:numref:`globalPatterDofManagerFig` shows the global pattern with the field-based ordering of unknowns.
Different colors mean different fields.
Red unkwnons are associated with displacement, yellow ones with pressure and blue ones with mass matrix.
Orange means the coupling among displacement and pressure, while green is the symmetric coupling among pressure and mass matrix.

.. _globalPatterDofManagerFig:
.. figure:: /coreComponents/linearAlgebra/docs/images/global.svg
   :align: center
   :width: 400
   :figclass: align-center

   Global pattern with field-based ordering.
   Red is associated with displacement unknowns, yellow with pressure ones and blue with those of mass matrix field.
   Orange means the coupling among displacement and pressure, while green is the symmetric coupling among pressure and mass matrix.

:numref:`permutedPatterDofManagerFig` shows the global pattern with the MPI rank-based ordering of unknowns.
In this case, just two processes are used.
Again, different colors indicate different ranks.

.. _permutedPatterDofManagerFig:
.. figure:: /coreComponents/linearAlgebra/docs/images/permutedGlobal.svg
   :align: center
   :width: 400
   :figclass: align-center

   Global pattern with MPI rank-based ordering.
   Red unkwnons are owned by rank 0 and green ones by rank 1.
   Blue indicates the coupling among the two processes.
