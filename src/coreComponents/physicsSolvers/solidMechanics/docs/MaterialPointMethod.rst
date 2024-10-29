.. _SolidMechanicsMPM:

#####################################
Material Point Method Solver
#####################################

List of Symbols
===================

.. math::
   i &\equiv \text{grid node index} \notag \\
   p &\equiv \text{particle index} \notag \\
   t &\equiv \text{time} \notag \\
   m &\equiv \text{mass} \notag \\
   v &\equiv \text{velocity} \notag \\
   V &\equiv \text{volume} \notag \\
   E &\equiv \text{Young's modulus} \notag \\
   \nu &\equiv \text{Poisson's ratio} \notag \\
   K &\equiv \text{bulk modulus} \notag \\
   G &\equiv \text{shear modulus} \notag \\
   \rho &\equiv \text{density} \notag \\
   \sigma &\equiv \text{stress} \notag \\
   P &\equiv \text{pressure} \notag \\
   c &\equiv \text{wavespeed} \notag \\
   \textbf{F} &\equiv \text{deformation gradient} \notag \\
   \textbf{L} \text{velocity gradient} \notag \\
   \mu &\equiv \text{friction coefficient} \notag \\

Introduction
==============
The `SolidMechanics_MPM` solver performs explicit time step simulations in which simulation data is carried by Lagrangian particles which are interpolated to a Eulerian background grid to solve linear momentum balance equation.

#. The material point method (MPM) tracks the material state at Lagrangian particles, the MPM avoids the advection errors that occur in Eulerian codes, making it possible to use constitutive models that have complex history dependence, as in many models of brittle damage.
#. Because the equation of mtion is solved on a fixed background grid, the MPM avoids the mesh entanglement that can occur in Lagrangian FEM, making the MPM suitbale for large-deformation solid mechanics.
#. The MPM has a natural no-slip no-interpenetration contact condition, and can robustly simulate frictional contact even with large deformation.
#. A MPM problem is initialized by defining a set of particles, which can be constructed directly from voxelized data (e.g. micro-computed tomography scans) without a costly mesh generation setup.

These benefits make the MPM appealing for mesoscale simulation of porous, heterogeneous, brittle, and granular materials.

The GEOSX-MPM solver implementation employs an uniformly spatially partitioned background grid algorithmically generated from user input paramters. Each partition corresponds to a parallel process in which adjacent partitions include an extra layer of ghost cells from their neighbors.
Problems are initilaized from an xml input which defines parameters of the background mesh (grid), constitutive models, boundary conditions, and output options. To facility easier job setup, a companion particle file writer python program can be used to define input parameters including the geometry and generate the xml input file and particle file required to run simulations.


MPM Algorithm
=========================

Particle-To-Grid Mapping
--------------------------
At the start of the time step (:math:`t = t^n`),  the mass and momentum are mapped from the partciles to a background grid. 
The mass (:math:'m_i') and velocity (:math:'v_i') for the :math:'i_^{th}' node of the background grid are:

.. math::
   m_{i}^{n}  = \sum_p S_{i,p} m_{p}^{n},

.. math::
   m_{i}^{n}  = \frac{1}{m_{i}^{n} } \sum_p S_{i,p} m_{p}^{n} v_{p}^{n},

where :math:`m_p` and :math:`v_p` are the mass and velocity of the :math:`p^{th}` particle, and :math:`S_{i,p}` is the average of the :math:`i^{th}` grid shape function :math:`S_i` over the :math:`p^{th}` particle domain. In our work, we use the convected particle domain interpolation (CPDI) formulation Sadeghriad et al. [2011] along with a domain scaling modification Homel et al. [2015] that allows for efficient parallelization and prevents nonphysical stretching or numerical fracture of failed particles.
The background grid is a structured, Cartesian mesh, with tri-linear shape functions, which simplifies evaluation of :math:`S_i`.

Internal and External Forces
-----------------------------

The internal and external force (neglecting surface tractions) at the :math:'i^{th}' grid node are:

.. math::
   f_{i}^{ext,(n)}  = \sum_p S_{i,p} b_p^{(n)} m_{p}^{(n)},

.. math::
   f_{i}^{int,(n)}  = -\sum_p \nabla S_{i,p} ( \textbf{\sigma}_p^{(n-1/2)} - q_{p}^{(n-1/2)} \textbf{I} ) V_p^{(n-1/2)},
   
where :math:'b_p' is the body force acting on the :math:'p^{th}' particle, :math:'V_p' is the particle volume, :math:'\textbf{\sigma}_p' is the Cauchy stress at the particle, and :math:'\bar{\nabla S_{i,p}}' is the average of the gradient of the :math:'i^{th}' grid shape function over the :math:'p^{th}' particle domain and :math:'q_p' is the artificial viscosity.

Trial Grid Velocity and Acceleration
------------------------------------

The trial (before contact) accelerations are computed on this background grid, from which the trial grid-node velocities are computed.

.. math::
   \bar{a}_i^{(n)} = \frac{f_i^{\text{int},(n)} + f_i^{\text{ext},(n)}}{m_i^{(n)}

Now we integrate acceleration at the grid to compute \bar{v}_i^{(n+1)}, the "trial" updated velocity, neglecting the effect of contact forces.

.. math::
   \bar{v}_i^{(n+1)} = v_i^{(n)} + \bar{a}_i^{(n)} \Delta t

Contact Forces
--------------------------
Contact forces (described below) modify the velocity field to compute a corrected :math:'v_i^{n+1}' from :math:'\bar{v}_i^{(n+1)}', and a corresponding chagen to the acceleration

.. math::
   v_i^{(n+1)} = \bar{v}_i^{(n+1)} + \frac{f_i^{\text{contact,(n+1)}}{m_{i}^{(n)}} \Delta t

.. math::
   a_i^{(n+1/2)} = \bar{a}_i^{(n)} + \frac{f_i^{\text{contact,(n+1)}}{m_{i}^{(n)}}

Cohesive Zones
--------------------------

In development (Will be updated later)

Grid Velocity Boundary Conditions
----------------------------------

Particle Kinematics
--------------------------

The updated grid velocity field is used to compute the updated velocity gradient at the particle, but here we use the half-step velocity

.. math::
   \textbf{L}_p^{(n+1/2)}=\sum_i ( v_{i}^{(n+1)-\frac{1}{2}a_{i}^{(n+1/2)}}) \tens \bar{\nabla S_{i,p}}

It is observed that a boundary instability arises when using the prescribed velocity gradient boundary condition with this half-step update, and the results were more stable when using :math:'\textbf{L}_p^{(n+1)} = \sum_i v_{i}^{(n+1) \tens \bar{\nabla S_{i,p}}'.
This option is the default when using either the prescribed deformation gradient table (fTable) to drive domain or boundary deformation.

The updated particle velocity gradient is used to update the particle deformation gradient,

.. math::
   \textbf{L}_p^{(n+1/2)} = \sum_i (v_{i}^{(n+1)} - \fracs{1}{2} a_{i}^{(n+1/2)} \Delta t) \tens \bar{\nabla S_{i,p}}

This integration is currently done using sub-stepping unless :code-block:'Fsubcycles="1"'. The updated deformation gradient is used to update the Jacobian, volume, and density:

.. math::
   \textbf{F}_p^{(n+1/2)} = \textbf{F}_p^{(n-1/2)} + (\textbf{L}_p^{(n+1/2)} \cdot \textbf{F}_p^{(n-1/2)}) \Delta t

Here, if using :code-block:'overlapCorrection="2"', we apply overlap correction to scale :math:'\textbf{F}_p^{(n+1/2)}' and :math:'J^{(n+1/2)}'. After that, we compute the updated volume and density:

.. math::
   V_{p}^{(n+1/2)} = V_{p}^{(0)} J_{p}^{(n+1/2)}

.. math::
   \rho_{p}^{(n+1/2)} = m_{p}/V_{p}^{(n+1/2)}

Artificial Viscosity
--------------------------

Internal Energy (Part 1)
--------------------------

Constitutive Update
--------------------------

Internal Energy (Part 2)
--------------------------

Grid-To-Particle Mapping
--------------------------

Fluid Implicit Particle (FLIP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The updated nodal velocities and accelerations are used to update the position and velocity of the particles. The default behavior is FLIP (:code-block:'updateMethod="FLIP"'):

.. math::
   x_{p}^{(n+1)} = x_{p}^{(n)} + \sum_{i} \bar{S}_{i,p}(v_i^{(n+1)-\frac{1}{2} a_i^{(n)} \Delta t) \Delta t

.. math::
   v_{p}^{(n+1)} = v_{p}^{(n)} + \sum_{i} \bar{S}_{i,p} a_i^{(n)} \Delta t

Particle-In-Cell (PIC)
^^^^^^^^^^^^^^^^^^^^^^^^^^

The more dissipation and stable option is PIC (:code-block:'updateMethod="PIC"'):

.. math::
   x_{p}^{(n+1)} = x_{p}^{(n)} + \sum_{i} \bar{S}_{i,p}(v_i^{(n+1)-\frac{1}{2} a_i^{(n)} \Delta t) \Delta t

.. math::
   v_{p}^{(n+1)} = v_{p}^{(n)} + \sum_{i} \bar{S}_{i,p} v_i^{(n+1)


XPIC
^^^^^^^^^^^^^^^^^^^^^^^^^^

Full Mass Matrix Material Point Method (FMPM)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Grid Update
--------------------------

Now, when using the prescribed deformation gradient table boundary condition the grid node locations and spcaings are updated, consistent with the prescribed velocity gradient.

.. math::
   x_{i}^{(n+1)} = x_{i}^{(n)} + \textbf{L}_\Gamma^{(n+1/2)} x_i^{(n)} \Delta T

.. math::
   [\Delta x, \Delta y, \Delta z]^{(n+1)} = [\Delta x, \Delta y, \Delta z]^{(n)} + \textbf{L}_\Gamma^{(n+1/2)} \cdot [\Delta x, \Delta y, \Delta z]^{(n)} \Delta t

Update Stable Time-Step
--------------------------

Each time-step iteration, the wavespeed is estimated for every particle by:

.. math:::
   c_p = \sqrt{\frac{K'_p^{(n+1/2) + \frac{4}{3}G'_p^{(n+1/2)}}{\rho_p^(n+1/2)}}

where :math:'K'' and math:'G'' are the effective bulk and shear moduli. The properties are termed effective because anisotropic materials have wavespeeds that vary with direction. In these instances, the direction with the maximum wavespeed is used.
The maximum value of :math:'c_p' in the domain is then used to compute the updated stable time-step:

.. math::
   ( \Delta t )^{n+1} = t^{n+2}-t^{n+1} = \text{CFL} \frac{min[\Delta x^{(n+1)},\Delta y^{(n+1)},\Delta z^{(n+1)}]}{max[c_p]+max[||v_p^{(n+1)||}]}

Fracture
=========================
A continuum-damage fracture treatment using the Field-Gradient Partitioning method is available, using either of the brittle damage laws described in Homel and Herbold [2017].
Additionally, a Weibull distribution of intiial strength ascribed to Voronoi cell groupings of particles can be defined to give more realistic brittle fracture response.
The Voronoi cell grouping of particle strength scales ensures correct nucleation and propagation of cracks.

.. code-block::

   damageFieldPartitioning="1"
   separabilityMinDamage="0.5"
   optimizeBinSort="0"
   binSizMultipllier="4"
   neighborRadius="1.0"

Grid Boundary Conditions
=========================

Boundary conditions for each face of the domain must be specified. The boundary condition types are:
| 0 = outflow
| 1 = symmetry
| 2 = boundary driven deformation (equivalent to a frictional rigid piston)

For type 2 boundary conditions, a deformation table (fTable) corresponding to the deformation of the background grid must be specified by a :math:'N"x4 matrix where :math:'N' is the number of rows. Each row begins with the time followed by the fractional deformation of the grid in x, y, and z directions.
The boundary velocity enforced is relative to the position with respect to the origin (:math:'x') given by :math:'v = L \cdot x' where :math:'L' is the velocity gradient calculated from the deformation table at a given timestep. Therefore, boundary velocities will be different for domains positioned differently with respect to the origin.
Note the solver variable prescribedBoundaryFTable must be set to 1 for the deformation table to take effect and the first row of the deformation table must always have deformations of 1 for each direction.

Prescribed transverse velocities can be ascribed to each boundary face, such as to rigidly affix the end of a cantilever, by specifiying :code-block:'enablePrescribedBoundaryTransverseVelocities=1' in the input file.

Additionally, periodic boundary conditions can be specified for any pair of opposing faces by an array of size 3 corresponding to X, Y, Z faces in that order. For example, the follling specifies periodic boundaries in the x direction by not y or z directions.

.. code-block::

   periodic="{1, 0, 0}"

Note that the corresponding boundary condition type for any periodic direction must be 0 (outflow) for particles to traverse the boundary. 
Periodic boundaries are useful for simulating RVEs where the domain would otherwise need to be large to negate boundary effects such as simulating the unit cell of a periodic architecture or microstructure.

A stress-driven boundary condition is implemented by a PID controller using the average domain stress to driven the domain deformation gradient. Care must be taken to choose PID gains that provide stable conditions. 
Likewise, some domains such as those that may fracture and have zero strength could become unstable.

Body Loads
=========================

By default there is no body force but a uniform one can be defined by specifying a 3 element array in the input file (e.g. :code-block:'bodyForce={b_x, b_y, b_z}').
The body forces for performing the generalized vortex problem for the method of manufactured solutions is enabled by specifiying :code-block:'generalizedVortexMMS="1"' in the input file.


Constitutive Models
=========================

Hyperelastic
--------------------------

A hyperelastic compressible neo-Hookean model is used to compute the elastic trial stress Fu and Ogden [2001].

.. math::
   \sigma^E = (\lambda \frac{ln J}{J} - \frac{\mu}{J}) \textbf{I} + \frac{\mu}{J} \textbf{F} \cdot \textbf{F}^T 

Hyperelastic MMS
--------------------------
A hyperelastic model of the form:


which is included for the method of manufactured solutions generalized vortex gradient problem.

ElasticIsotropic
--------------------------
Hypoelastic elastic isotropic constitutive model

Perfectly Plastic
--------------------------
Hypoelastic elastic isotropic constitutive model defined by the material density, any two elastic constants (e.g. Young's modulus and Poisson's ratio or bulk and shear moduli ) and a yield strength.


VonMisesJ
--------------------------
A model identical to perfectly plastic except that the pressure is computed exactly from the instantaneous Particle jacobian given by :math:'p = -K log{J_{p}}' rather than a hypoelastic update. This was implemetned for verification testing of the overlap correction code (which scales J) but can be used more generally.
A tpyical paramterization for the model is as follows, where stress has the units of GPa and density is mg/mm3.

.. code-block::

   <VonMisesJ
      name="VonMisesJ"
      defaultDensity="8.96"
      defaultBulkModulus="130"
      defaultShearModulus="44"
      defaultYieldStrength="0.322"/>


StrainHardeningPolymer
-------------------------

CeramicDamage
-------------------------

Inputs are:

| :math:'\rho' = material density
| :math:'K' = bulk modulus
| :math:'G' = shear modulus
| :math:'\sigma_t' = reference tensile strength
| :math:'\sigma_c' = reference compressive strength
| :math:'\sigma_{max}' = maximum compressive strength
| :math:'\mu' = shear frictional slope for fully damaged granular material
| = crack speed

.. code-block::

   <CeramicDamage
	 name="sand"
	 defaultDensity="2.648"
	 defaultBulkModulus="36.3"
	 defaultShearModulus="26.0"
	 tensileStrength="0.449"
	 compressiveStrength="2.27"
	 maximumStrength="5.0"
	 crackSpeed="1.8"
	 thirdInvariantDependence="1"/>


Graphite
--------------------------

Transversely isotropic constitutive model defined by a basal plane material direction

Geomechanics
-------------------------


Weibull Variability
=========================
When using a smeared (continuum) damage model to simulate brittle failure, it is necessary to regularize the solution to avoid mesh dependence, for example by introducing spatial heterogeneity to decsribe the aleatory uncertainty in the material strength (c.f. Strack et al. [2014]).
A classical approach is to define the particle strength, :math:'\sigma_f', to be both spatially variable and size-dependent according to a Weibull distribution (Weibull [1951], Jayatilaka and Trustrum [1977]).

.. math::
   \frac{\sigma_f}{\bar{\sigma}_f} = (\frac{\bar{V} ln R }{V ln 1/2})^{1/m}

where :math:'V' is the particle volume, :math:'\bar{\sigma}_f' and :math:'\bar{V}' ar the reference failure stress and corresponding reference volume,  :math:'m' is the Weibull modulus, and :math:'R' is a uniformly distributed random number, :math:'0<R<1'.
The 2-D computational examples in this work are all plane-strain, but it is necessary to specify a material thickness to define the particle volume. While there is some evidence that a Weibull distribution may not be physically representative of flaw distributions in many engineering materials (c.f. Tonge and Ramesh [2016]), the need for physically motivated spatial heterogeneity is clear.

When using Weibull Variability you need to specify a reference volume and Weibull modulus that will be sued for all materials that have Weibull scaled strength.
The particle file write geometry object voronoiWeibullBoxWrapper can be used to generate strength scaling for particles during preprocessing of the particle file.

Particle and Grid fields
=========================
.. list-table:: Particle fields

   * - Field Name
     - Dimensions
     - Description
   * - particleMass
     - {num particles}
     - mass of each particle
   * - particleVolume
     - {num particles}
     - volume of each particle
   
.. list-table:: Grid fields

   * - Field Name
     - Dimensions
     - Description
   * - gridMass
     - {num particles, num velocity fields}
     - mass of each field at the grid node


Solver Events
=============

The MPM solver supports its own events for various functionallities. Each event 

Material Swap
-----------------------------------

Performs a material swap essentially copying particles from one particle region to another. Input file must include an empty particle region with the desired material model to copy particles to.

.. code-block::
   
   <ElasticIsotropic
	 name="mat1"
	 defaultDensity="2.7"
	 defaultBulkModulus="6.7"
	 defaultShearModulus="2.6"/>

   <ElasticIsotropic
	 name="mat2"
	 defaultDensity="2.7"
	 defaultBulkModulus="67"
	 defaultShearModulus="26"/>

.. code-block::

   <MPMEvents>
      <MaterialSwap 
       time="1.0"
       interval="0.1"
       source="ParticleRegion1"
       destination="ParticleRegion2"/>
   </MPMEvents>

Anneal
-----------------------------------

.. code-block::
  
   <MPMEvents>
      <Anneal 
       time=1.0
       interval=1.0
       source="all" />
   </MPMEvents>

Heal
-----------------------------------

Crystal Heal
-----------------------------------

Insert Periodic Contact Surfaces
-----------------------------------

Machine Sample
-----------------------------------

Friction Coefficient Swap
-----------------------------------

Body Force Update
-----------------------------------

Deformation Update
-----------------------------------

Cohesive Zone Reference
-----------------------------------

Solver Specific Parameters
==========================

This solver has many solver specific variables. The following list describes each and their use.

.. list-table:: Solution Control
   
   * - type
     - input
     - default
     - description
   * - bool
     - solverProfiling
     - 0
     - output solver specific timing information to console
   * - array1d< string >
     - plottableFields
     - all fields
     - a string array of particle and grid fields to write to vtk output files
   * - string
     - timeIntegrationOption
     - ExplicitDynamic
     - Time integration option (QausiStatic, ImplicitDynamic, ExplicitDyanmic)
   * - real 
     - cflFactor
     - 0.25
     - Courant number 
   * - real
     - initialDt
     - 1e-16
     - Initial time step
   * - string
     - updateMethod
     - FLIP
     - Grid-to-particle update method (FLIP, PIC, XPIC, FMPM)
   * - int
     - updateOrder
     - 2
     - Update order for XPIC and FMPM update methods (ignored by FLIP and PIC)
   * - array1d< real >
     - bodyForce
     - {0.0, 0.0, 0.0}
     - Uniform body force with components :math:'b_x', :math:'b_y', :math:'b_z'
   * - bool
     - boxAverageHistory
     - 0
     - Enables output of box average values
   * - bool
     - reactionHistory
     - 0
     - Enables output of global boundary reactions
   * - bool
     - planeStrain
     - 0
     - Enables plane strain calculations for 2-D simulations
   * - bool
     - plotUnscaledParticles
     - 0
     - Determines whether particle domains are scaled when plotting (CPDI domain scaling)
   * - bool
     - cpdiDomainScaling
     - 0 
     - Enables CPDI domain scaling
   * - bool
     - debugFlag
     - 0
     - Enables writing of solver function calls to console (useful for debugging in release compilation)
   * - bool
     - generalizedVortexMMS
     - 0
     - Enables the body forces for performing the generalize vortex problem for the method of manufactured solutions
   * - bool
     - needsNeighborList
     - 0
     - Enables neighbor list generation (will be automatically set for any SPH operations)

.. list-table:: Field Gradient Partitioning

   * - type
     - input
     - default
     - description
   * - bool
     - damageFieldPartitioning
     - 0
     - set to 1 for fracture or self-contact
   * - real
     - separabilityMinDamage
     - 0.5
     - mininum average damage for separability
   * - bool
     - treatFullyDamagedAsSingleField
     - 0
     - will treat regions of full damage as a single field (recommended)


.. list-table:: Artificial Viscosity

   * - type
     - input
     - default
     - description
   * - bool
     - useArtificialViscosity
     - 0
     - Enables use of artificial viscosity
   * - real
     - artificialViscosityQ0
     - 0.0
     - Linear term, 1.5 is a common value
   * - real
     - artificialViscosityQ1
     - 0.0
     - Quadratic term, 0.5 is a common value

.. list-table:: Boundary Conditions

   * - type
     - input
     - default
     - description
   * - bool
     - prescribedBcTable
     - 0
     - Determines if a boundary condition table is used
   * - array2d< real >
     - bcTable
     - empty
     - A table where column 1 defines the simulation time and columns 2-7 define the boundary condition type for -x, +x, -y, +y, -z, and +z faces respectively. Each row denotes a different time and set of boundary conditions.
   * - array1d< real >
     - boundaryConditionTypes
     - {0,0,0,0,0,0}
     - A 6 element vector defining the boundary condition types for faces -x, +x, -y, +y, -z, and +z faces
   * - bool
     - prescribedBoundaryFTable
     - 0
     - Enables boundary driven deformation using interpolated fTable for faces with boundary condition type 2
   * - bool
     - prescribedFtable
     - 0
     - Enables superimposed velocity gradient for quasistatic loading under triply periodic boundary conditions
   * - string
     - fTableInterpType
     - Linear
     - Determines the interpolation method for the fTable (Linear, Cosine, Smoothstep )
   * - array2d< real >
     - fTable
     - {{0.0, 1.0, 1.0, 1.0}}
     - A table where column 1 defines the simulation time and columns 2-4 define the deformation value of the domain in the x, y, and z directions respectively. Each row denotes a different time and deformation values.
   * - bool
     - enablePrescribedBoundaryTransverseVelocities
     - 0
     - Enables prescribed transverse velocities on boundary faces
   * - array2d< real >
     - prescribedBoundaryTransverseVelocities
     - {{0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}}
     - A 6x2 table that defines the transverse velocity bondary conditions for faces -x, +x, -y, +y, -z, and +z respectively.
   * - array1d< bool >
     - stressControl
     - {0, 0, 0}
     - A vector which enables PID stress control for x, y, and z directions.
   * - real
     - stressControlKp
     - 0.0
     - P gain of PID controller for stress control
   * - real
     - stressControlKi
     - 0.0
     - I gain of PID controller for stress control
   * - real
     - stressControlKd
     - 0.0
     - D gain of PID controller for stress control

.. list-table:: Cohesive Zones

   * - type
     - input
     - default
     - description

.. list-table:: Events

   * - type
     - input
     - default
     - description
   * - bool
     - useEvents
     - 0 
     - Enables use of MPM solver events

.. list-table:: Contact

   * - type
     - input
     - default
     - description
   * - real
     - explicitSurfaceNormalInfluence
     - 0.0
     - Value determines weighting of explicit particle surface normals from implicit grid surface normals for contact calculations
   * - bool
     - useSurfacePositionForContact
     - 0
     - Enables use of explicit particle surface positions for more accurate contact gap calculations (may be unstable for severe deformation)
   * - real
     - frictionCoefficient
     - 0.0
     - Global friction coefficient value, overriden if friction coefficient table is specified
   * - array2d< real >
     - frictionCoefficientTable
     - 
     - 

Examples
=========================

An example of a valid XML block is given here:

.. literalinclude:: ../../../../../inputFiles/solidMechanics/sedov_finiteStrain_smoke.xml
  :language: xml
  :start-after: <!-- SPHINX_SOLID_MECHANICS_SOLVER -->
  :end-before: <!-- SPHINX_SOLID_MECHANICS_SOLVER_END -->


How to Cite
==========================

References
=========================