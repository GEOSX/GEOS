.. _ExampleMCCWellbore:


####################################################
Modified Cam-Clay Model for Wellbore Problems
####################################################


**Context**

In this benchmark example, the Modified Cam-Clay model (see :ref:`ModifiedCamClay`) is applied to solve for elastoplastic deformation within the vicinity of a vertical wellbore. For the presented example, an analytical solution is employed to verify the accuracy of the numerical results. The resulting model can be used as a base for more complex analysis (e.g., wellbore drilling, fluid injection and storage scenarios). 


**Input file**

Everything required is contained within two GEOS input files located at:

.. code-block:: console

  inputFiles/solidMechanics/ModifiedCamClayWellbore_base.xml

.. code-block:: console

  inputFiles/solidMechanics/ModifiedCamClayWellbore_benchmark.xml

------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

We simulate a drained wellbore problem subjected to isotropic horizontal stress (:math:`\sigma_h`) and vertical stress (:math:`\sigma_v`), as shown below. By increasing the wellbore supporting pressure (:math:`P_w`), the wellbore expands, and the formation rock experiences elastoplastic deformation. A plastic zone develops in the near wellbore region.


.. _problemSketchMccWellFig:
.. figure:: MCC_wellSketch.png
   :align: center
   :width: 500
   :figclass: align-center

   Sketch of the wellbore problem 


To simulate this phenomenon, the Modified Cam-Clay model is used in this example. Displacement and stress fields around the wellbore are numerically calculated. These numerical predictions are then compared with the corresponding analytical solutions `(Chen and Abousleiman, 2013)  <https://www.icevirtuallibrary.com/doi/10.1680/geot.11.P.088>`__ from the literature. 



In this example, we focus our attention on the ``Mesh`` tags,
the ``Constitutive`` tags, and the ``FieldSpecifications`` tags.

------------------------------------------------------------------
Mesh
------------------------------------------------------------------

Following figure shows the generated mesh that is used for solving this wellbore problem.


.. _problemSketchMccMeshFig:
.. figure:: MCC_wellMesh.png
   :align: center
   :width: 500
   :figclass: align-center

   Generated mesh for a vertical wellbore problem


Let us take a closer look at the geometry of this wellbore problem.
We use the internal wellbore mesh generator ``InternalWellbore`` to create a rock domain
(:math:`10\, m \, \times 5 \,  m \, \times 2 \, m`), with a wellbore of
initial radius equal to :math:`0.1` m.
Coordinates of ``trajectory`` defines the wellbore trajectory, which represents a vertical well in this example. 
By turning on ``autoSpaceRadialElems="{ 1 }"``, the internal mesh generator automatically sets number and spacing of elements in the radial direction, which overrides the values of ``nr``. 
With ``useCartesianOuterBoundary="0"``, a Cartesian aligned outer boundary on the outer block is enforced.
In this way, a structured three-dimensional mesh is created with 50 x 40 x 2 elements in the radial, tangential and z directions, respectively. All the elements are eight-node hexahedral elements (``C3D8``) and refinement is performed
to conform with the wellbore geometry. This mesh is defined as a cell block with the name
``cb1``.


.. literalinclude:: ../../../../../../../inputFiles/solidMechanics/ModifiedCamClayWellbore_benchmark.xml
    :language: xml
    :start-after: <!-- SPHINX_WELLBORE_MESH -->
    :end-before: <!-- SPHINX_WELLBORE_MESH_END -->


------------------------
Solid mechanics solver
------------------------

For the drained wellbore problem, the pore pressure variation is omitted. Therefore, we just need to define a solid mechanics solver, which is called ``mechanicsSolver``. 
This solid mechanics solver (see :ref:`SolidMechanicsLagrangianFEM`) is based on the Lagrangian finite element formulation. 
The problem is run as ``QuasiStatic`` without considering inertial effects. 
The computational domain is discretized by ``FE1``, which is defined in the ``NumericalMethods`` section. 
The material is named as ``rock``, whose mechanical properties are specified in the ``Constitutive`` section.

.. literalinclude:: ../../../../../../../inputFiles/solidMechanics/ModifiedCamClayWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_WELLBORE_SOLVER -->
  :end-before: <!-- SPHINX_WELLBORE_SOLVER_END -->


------------------------------
Constitutive laws
------------------------------

For this drained wellbore problem, we simulate the elastoplastic deformation caused by wellbore expansion.
A homogeneous domain with one solid material is assumed, whose mechanical properties are specified in the ``Constitutive`` section: 

.. literalinclude:: ../../../../../../../inputFiles/solidMechanics/ModifiedCamClayWellbore_base.xml
    :language: xml
    :start-after: <!-- SPHINX_WELLBORE_MATERIAL -->
    :end-before: <!-- SPHINX_WELLBORE_MATERIAL_END -->


Recall that in the ``SolidMechanics_LagrangianFEM`` section, 
``rock`` is designated as the material in the computational domain. 
Here, Modified Cam-Clay ``ModifiedCamClay`` is used to simulate the elastoplastic behavior of ``rock``.

The following material parameters should be defined properly to reproduce the analytical example:

.. include:: /docs/sphinx/datastructure/ModifiedCamClay.rst


The constitutive parameters such as the density, the bulk modulus, and the shear modulus are specified in the International System of Units.


-----------------------------------------------------------
Initial and boundary conditions
-----------------------------------------------------------

The next step is to specify fields, including:

  - The initial value (the in-situ stresses and traction at the wellbore wall have to be initialized)
  - The boundary conditions (the reduction of wellbore pressure and constraints of the outer boundaries have to be set)

In this tutorial, we need to specify isotropic horizontal stress (:math:`\sigma_h` = -100 kPa) and vertical stress (:math:`\sigma_v` = -160 kPa). 
To reach equilibrium, a compressive traction :math:`P_w` = -100 kPa is instantaneously applied at the wellbore wall ``rneg`` at time :math:`t` = 0 s, which will then be gradually increased to a higher value (-300 kPa) to let wellbore expand.
The remaining parts of the outer boundaries are subjected to roller constraints.  
These boundary conditions are set up through the ``FieldSpecifications`` section.


.. literalinclude:: ../../../../../../../inputFiles/solidMechanics/ModifiedCamClayWellbore_base.xml
    :language: xml
    :start-after: <!-- SPHINX_WELLBORE_BC -->
    :end-before: <!-- SPHINX_WELLBORE_BC_END -->


With ``tractionType="normal"``, traction is applied to the wellbore wall ``rneg`` as a pressure specified from the product of scale ``scale="-1.0e5"`` and the outward face normal. 
A table function ``timeFunction`` is used to define the time-dependent traction ``ExternalLoad``. 
The ``coordinates`` and ``values`` form a time-magnitude
pair for the loading time history. In this case, the loading magnitude increases linearly as the time evolves. 

.. literalinclude:: ../../../../../../../inputFiles/solidMechanics/ModifiedCamClayWellbore_base.xml
    :language: xml
    :start-after: <!-- SPHINX_WELLBORE_TABLE -->
    :end-before: <!-- SPHINX_WELLBORE_TABLE_END -->

You may note :

 - All initial value fields must have ``initialCondition`` field set to ``1``;
 - The ``setName`` field points to the previously defined set to apply the fields;
 - ``nodeManager`` and ``faceManager`` in the ``objectPath`` indicate that the boundary conditions are applied to the element nodes and faces, respectively;
 - ``fieldName`` is the name of the field registered in GEOS;
 - Component ``0``, ``1``, and ``2`` refer to the x, y, and z direction, respectively;
 - And the non-zero values given by ``scale`` indicate the magnitude of the loading; 
 - Some shorthand, such as ``xneg`` and ``xpos``, are used as the locations where the boundary conditions are applied in the computational domain. For instance, ``xneg`` means the face of the computational domain located at the left-most extent in the x-axis, while ``xpos`` refers to the face located at the right-most extent in the x-axis. Similar shorthands include ``ypos``, ``yneg``, ``zpos``, and ``zneg``;
 - The mud pressure loading has a negative value due to the negative sign convention for compressive stress in GEOS. 

 
The parameters used in the simulation are summarized in the following table.

+------------------+---------------------------+------------------+---------------+
| Symbol           | Parameter                 | Units            | Value         |
+==================+===========================+==================+===============+
| :math:`P_r`      | Reference Pressure        | [kPa]            | 120           |
+------------------+---------------------------+------------------+---------------+
| :math:`G`        | Shear Modulus             | [kPa]            | 4302          |
+------------------+---------------------------+------------------+---------------+
| :math:`P_c`      | PreConsolidation Pressure | [kPa]            | 169           |
+------------------+---------------------------+------------------+---------------+
| :math:`M`        | Slope of CSL              | [-]              | 1.2           |
+------------------+---------------------------+------------------+---------------+
| :math:`c_c`      | Virgin Compression Index  | [-]              | 0.072676      |
+------------------+---------------------------+------------------+---------------+
| :math:`c_r`      | Recompression Index       | [-]              | 0.014535      |
+------------------+---------------------------+------------------+---------------+
| :math:`\sigma_h` | Horizontal Stress         | [kPa]            | -100          |
+------------------+---------------------------+------------------+---------------+
| :math:`\sigma_v` | Vertical Stress           | [kPa]            | -160          |
+------------------+---------------------------+------------------+---------------+
| :math:`a_0`      | Initial Well Radius       | [m]              | 0.1           |
+------------------+---------------------------+------------------+---------------+
| :math:`P_w`      | Mud Pressure              | [kPa]            | -300          |
+------------------+---------------------------+------------------+---------------+


---------------------------------
Inspecting results
---------------------------------

In the above example, we requested silo-format output files. We can therefore import these into VisIt and use python scripts to visualize the outcome. The following figure shows the distribution of :math:`\sigma_{\theta\theta}` in the near wellbore region.

.. _problemVerificationMccFig1:
.. figure:: MCC_Stress.png
   :align: center
   :width: 1000
   :figclass: align-center

   Simulation result of :math:`\sigma_{\theta\theta}`



The figure below shows the comparisons between the numerical predictions (marks) and the corresponding analytical solutions (solid curves) with respect to the distributions of normal stress components, stress path, the supporting wellbore pressure and wellbore size. It is evident that the predictions well match the analytical results.  



.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/mccWellbore/mccWellboreFigure.py



------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.






