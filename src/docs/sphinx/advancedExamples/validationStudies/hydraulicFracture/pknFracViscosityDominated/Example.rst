.. _pknViscosityDominated:


###########################################################
Viscosity-Storage-Dominated PKN Hydraulic Fracture
###########################################################


**Context**

In this example, we simulate the propagation of a Perkins–Kern–Nordgren (PKN) fracture in a viscosity-storage-dominated regime, a classic benchmark in hydraulic fracturing. The developed planar fracture displays an elliptical vertical cross-section. Unlike KGD and penny-shaped fractures, the growth height of a PKN fracture is constrained by mechanical barriers (such as bedding layers, sedimentary laminations, or weak interfaces), thus promoting lateral propagation. This problem is solved using the hydrofracture solver in GEOS to obtain the temporal evolutions of the fracture characteristics (length, aperture, and pressure). We validate these simulated values against existing analytical solutions `(Kovalyshen and Detournay, 2010;  <https://link.springer.com/article/10.1007/s11242-009-9403-4>`__ `Economides and Nolte, 2000)  <https://books.google.com/books/about/Reservoir_Stimulation.html?id=rDlQAQAAIAAJ>`__. 


**Input file**

This example uses no external input files. Everything we need is contained within two GEOS input files:

.. code-block:: console

  inputFiles/hydraulicFracturing/pknViscosityDominated_base.xml

.. code-block:: console

  inputFiles/hydraulicFracturing/pknViscosityDominated_benchmark.xml


Python scripts for post-processing and visualizing the simulation results are also prepared:

.. code-block:: console

  inputFiles/hydraulicFracturing/scripts/hydrofractureQueries.py

.. code-block:: console

  inputFiles/hydraulicFracturing/scripts/hydrofractureFigure.py


------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

In this example, a hydraulic fracture initiates and propagates from the center of a 20m-thick layer. This layer is homogeneous and bounded by neighboring upper and lower layers. For viscosity-dominated fractures, more energy is necessary to move the fracturing fluid than to split the intact rock. If fluid leak-off is neglected, storage-dominated propagation occurs with most of the injected fluid confined within the open surfaces. To meet the requirements of the viscosity-storage-dominated assumptions, impermeable domain (no fluid leak-off), incompressible fluid with constant viscosity (:math:`1.0 cp`) and ultra-low rock toughness (:math:`0.1 MPa { \sqrt{m} }`) are chosen in the GEOS simulation. With these parameters, the fracture stays within the target layer; it extends horizontally and meets the conditions of the PKN fracture in a viscosity-storage-dominated regime.

We assume that the fluid injected in the fracture follows the lubrication equation resulting from mass conservation and Poiseuille's law. The fracture propagates by creating new surfaces if the stress intensity factor exceeds the local rock toughness :math:`K_{IC}`. As the geometry of the PKN fracture exhibits symmetry, the simulation is reduced to a quarter-scale. For verification purposes, a plane strain deformation is considered in the numerical model.

We set up and solve a hydraulic fracture model to obtain the temporal solutions of the fracture half length :math:`l`, the net pressure :math:`p_0` and the fracture aperture :math:`w_0` at the fluid inlet for the PKN fracture propagating in viscosity-storage-dominated regime. `Kovalyshen and Detournay (2010)  <https://link.springer.com/article/10.1007/s11242-009-9403-4>`__ and `Economides and Nolte (2000)  <https://books.google.com/books/about/Reservoir_Stimulation.html?id=rDlQAQAAIAAJ>`__ derived the analytical solutions for this classic hydraulic fracture problem, used here to verify the results of the GEOS simulations:

.. math:: l(t) = 0.3817 (\frac{ E_p Q_0^3 t^4 }{ \mu h^4 })^{ 1/5 }

.. math:: w_0(t) = 3 (\frac{ \mu Q_0 l }{ E_p })^{ 1/4 } 

.. math:: p_0(t) = (\frac{ 16 \mu Q_0 E_p^3 l }{ \pi h^4 })^{ 1/4 } 

where the plane modulus :math:`E_p` is related to Young's modulus :math:`E` and Poisson's ratio :math:`\nu`:

.. math:: E_p = \frac{ E }{ 1-\nu^2 }



For this example, we focus on the ``Mesh``,
the ``Constitutive``, and the ``FieldSpecifications`` tags.

------------------------------------------------------------------
Mesh
------------------------------------------------------------------

The following figure shows the mesh used in this problem.


.. _problemSketchpknViscosityDominatedFig:
.. figure:: mesh.png
   :align: center
   :width: 500
   :figclass: align-center

   Generated mesh

We use the internal mesh generator to create a computational domain
(:math:`400\, m \, \times 400 \,  m \, \times 800 \, m`), as parametrized in the ``InternalMesh`` XML tag. 
The structured mesh contains 105 x 105 x 60 eight-node brick elements in the x, y, and z directions respectively. 
Such eight-node hexahedral elements are defined as ``C3D8`` elementTypes, and their collection forms a mesh
with one group of cell blocks named here ``cb1``. Local refinement is performed for the elements in the vicinity of the fracture plane. 


.. literalinclude:: ../../../../../../../inputFiles/hydraulicFracturing/pknViscosityDominated_benchmark.xml
    :language: xml
    :start-after: <!-- SPHINX_MESH -->
    :end-before: <!-- SPHINX_MESH_END -->


The fracture plane is defined by a nodeset occupying a small region within the computational domain, where the fracture tends to open and propagate upon fluid injection:


.. literalinclude:: ../../../../../../../inputFiles/hydraulicFracturing/pknViscosityDominated_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_FRACPLANE -->
  :end-before: <!-- SPHINX_FRACPLANE_END -->


------------------------
Solid mechanics solver
------------------------

GEOS is a multi-physics platform. Different combinations of
physics solvers available in the code can be applied
in different regions of the domain and be functional at different stages of the simulation.
The ``Solvers`` tag in the XML file is used to list and parameterize these solvers.

Three elementary solvers are combined in the solver ``hydrofracture`` to model the coupling between fluid flow within the fracture, rock deformation, fracture deformation and propagation:


.. literalinclude:: ../../../../../../../inputFiles/hydraulicFracturing/pknViscosityDominated_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACSOLVER -->
  :end-before: <!-- SPHINX_HYDROFRACSOLVER_END -->


- Rock and fracture deformations are modeled by the solid mechanics solver ``SolidMechanicsLagrangianSSLE``. In this solver, we define ``targetRegions`` that includes both the continuum region and the fracture region. The name of the contact constitutive behavior is specified in this solver by the ``contactRelationName``.


.. literalinclude:: ../../../../../../../inputFiles/hydraulicFracturing/pknViscosityDominated_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_MECHANICALSOLVER -->
  :end-before: <!-- SPHINX_MECHANICALSOLVER_END -->


- The single-phase fluid flow inside the fracture is solved by the finite volume method in the solver ``SinglePhaseFVM``.


.. literalinclude:: ../../../../../../../inputFiles/hydraulicFracturing/pknViscosityDominated_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_SINGLEPHASEFVM -->
  :end-before: <!-- SPHINX_SINGLEPHASEFVM_END -->


- The solver ``SurfaceGenerator`` defines the fracture region and rock toughness ``rockToughness="0.1e6"``. With ``nodeBasedSIF="1"``, a node-based Stress Intensity Factor (SIF) calculation is chosen for the fracture propagation criterion. 


.. literalinclude:: ../../../../../../../inputFiles/hydraulicFracturing/pknViscosityDominated_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_SURFACEGENERATOR -->
  :end-before: <!-- SPHINX_SURFACEGENERATOR_END -->


------------------------------
Constitutive laws
------------------------------

For this problem, a homogeneous and isotropic domain with one solid material is assumed. Its mechanical properties and associated fluid rheology are specified in the ``Constitutive`` section. 
``ElasticIsotropic`` model is used to describe the mechanical behavior of ``rock`` when subjected to fluid injection.
The single-phase fluid model ``CompressibleSinglePhaseFluid`` is selected to simulate the response of ``water`` upon fracture propagation.


.. literalinclude:: ../../../../../../../inputFiles/hydraulicFracturing/pknViscosityDominated_base.xml
    :language: xml
    :start-after: <!-- SPHINX_MATERIAL -->
    :end-before: <!-- SPHINX_MATERIAL_END -->


All constitutive parameters such as density, viscosity, bulk modulus, and shear modulus are specified in the International System of Units.


------------------------------
Time history function
------------------------------

In the ``Tasks`` section, ``PackCollection`` tasks are defined to collect time history information from fields. 
Either the entire field or specified named sets of indices in the field can be collected.
In this example, ``pressureCollection``, ``apertureCollection``, ``hydraulicApertureCollection`` and ``areaCollection`` are specified to output the time history of fracture characterisctics (pressure, width and area). 
``objectPath="ElementRegions/Fracture/FractureSubRegion"`` indicates that these ``PackCollection`` tasks are applied to the fracure element subregion.

.. literalinclude:: ../../../../../../../inputFiles/hydraulicFracturing/pknViscosityDominated_base.xml
    :language: xml
    :start-after: <!-- SPHINX_TASKS -->
    :end-before: <!-- SPHINX_TASKS_END -->

These tasks are triggered using the ``Event`` manager with a ``PeriodicEvent`` defined for the recurring tasks. 
GEOS writes one file named after the string defined in the ``filename`` keyword and formatted as a HDF5 file (``pknViscosityDominated_output.hdf5``). This TimeHistory file contains the collected time history information from specified time history collector.
This file includes datasets for the simulation time, fluid pressure, element aperture, hydraulic aperture and element area for the propagating hydraulic fracture.
A Python script is prepared to read and query any specified subset of the time history data for verification and visualization. 


-----------------------------------------------------------
Initial and boundary conditions
-----------------------------------------------------------

The next step is to specify:

  - The initial values: the ``waterDensity``, ``separableFace`` and the ``ruptureState`` of the propagating fracture have to be initialized,
  - The boundary conditions: fluid injection rates and the constraints of the outer boundaries have to be set.

In this example, a mass injection rate ``SourceFlux`` (``scale="-6.625"``) is applied at the surfaces of the initial fracture. Only one fourth of the total injection rate is used because only a quarter of the fracture is modeled (the problem is symmetric). The value given for ``scale`` is :math:`Q_0 \rho_f/4` (not :math:`Q_0 /4`). 
All the outer boundaries are subject to roller constraints. 
These boundary conditions are set through the ``FieldSpecifications`` section.


.. literalinclude:: ../../../../../../../inputFiles/hydraulicFracturing/pknViscosityDominated_base.xml
    :language: xml
    :start-after: <!-- SPHINX_BC -->
    :end-before: <!-- SPHINX_BC_END -->


The parameters used in the simulation are summarized in the following table.

+------------------+-------------------------+--------------------+--------------------+
| Symbol           | Parameter               | Unit               | Value              |
+==================+=========================+====================+====================+
| :math:`K`        | Bulk Modulus            | [GPa]              | 20.0               |
+------------------+-------------------------+--------------------+--------------------+
| :math:`G`        | Shear Modulus           | [GPa]              | 12.0               |
+------------------+-------------------------+--------------------+--------------------+
| :math:`K_{IC}`   | Rock Toughness          | [MPa.m\ :sup:`1/2`]| 0.1                |
+------------------+-------------------------+--------------------+--------------------+
| :math:`\mu`      | Fluid Viscosity         | [Pa.s]             | 1.0x10\ :sup:`-3`  |
+------------------+-------------------------+--------------------+--------------------+
| :math:`Q_0`      | Injection Rate          | [m\ :sup:`3`/s]    | 0.0265             |
+------------------+-------------------------+--------------------+--------------------+
| :math:`t_{inj}`  | Injection Time          | [s]                | 200                |
+------------------+-------------------------+--------------------+--------------------+
| :math:`h_f`      | Fracture Height         | [m]                | 20                 |
+------------------+-------------------------+--------------------+--------------------+



---------------------------------
Inspecting results
---------------------------------

The following figure shows the distribution of :math:`\sigma_{zz}` at :math:`t=200 s` within the computational domain..

.. _problemVerificationpknViscosityDominatedFig:
.. figure:: szz.png
   :align: center
   :width: 500
   :figclass: align-center

   Simulation result of :math:`\sigma_{zz}` at :math:`t=200 s`

First, by running the query script 

.. code-block:: console

   python ./hydrofractureQueries.py pknViscosityDominated

the HDF5 output is postprocessed and temporal evolution of fracture characterisctics (fluid pressure and fracture width at fluid inlet and fracure half length) are saved into a txt file ``model-results.txt``, which can be used for verification and visualization:

.. code-block:: console
		
  [['      time', '  pressure', '  aperture', '    length']]
           2 1.413e+06 0.0006093       5.6
           4 1.174e+06 0.0007132       8.4
           6 1.077e+06 0.0007849      10.8
           8 1.044e+06 0.0008482      12.8
          10 1.047e+06 0.0009098      14.8

Note: GEOS python tools ``geosx_xml_tools`` should be installed to run the query script (See `Python Tools Setup <https://geosx-geosx.readthedocs-hosted.com/projects/geosx-geospythonpackages/en/latest/>`_ for details). 
 
Next, figure below shows the comparisons between the results from GEOS simulations (markers) and the corresponding
analytical solutions (curves) for the example with viscosity-storage dominated assumptions, which is generated using the visualization script:

.. code-block:: console

   python ./pknViscosityDominatedFigure.py


The evolution in time of the fracture half-length, the near-wellbore fracture aperture, and the fluid pressure all correlate well with the analytical
solutions.  

.. plot:: docs/sphinx/advancedExamples/validationStudies/hydraulicFracture/pknFracViscosityDominated/pknViscosityDominatedFigure.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
