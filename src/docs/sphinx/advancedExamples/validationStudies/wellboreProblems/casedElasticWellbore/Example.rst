.. _AdvancedExampleCasedElasticWellbore:


####################################################
Cased Elastic Wellbore Problem
####################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example uses the solid mechanics solver to handle a cased wellbore problem subjected to a pressure test. The completed wellbore is composed of a steel casing, a cement sheath and rock formation. Isotropic linear elastic behavior is assumed for all the three materials. No separation is allowed for the casing-cement and cement-rock contact interfaces.

Analytical results of the radial and hoop stresses, :math:`\sigma_{rr}` and :math:`\sigma_{\theta\theta}`, in casing, cement sheath and rock are expressed as `(Hervé and Zaoui, 1995) <https://link.springer.com/chapter/10.1007%2F978-94-015-8494-4_55>`__ :

.. math::
   \sigma_{rr} = (2\lambda + 2G)A - \frac{2GB}{r^2}

.. math::
   \sigma_{\theta\theta} = (2\lambda + 2G)A + \frac{2GB}{r^2}

where :math:`\lambda` and :math:`G` are the Lamé moduli, :math:`r` is the radial coordinate, :math:`A` and :math:`B` are piecewise constants that are obtained by solving the boundary and interface conditions, as detailed in the post-processing script.


**Input file**

This benchmark example uses no external input files and everything required is
contained within two GEOS xml files that are located at:

.. code-block:: console

  inputFiles/wellbore/CasedElasticWellbore_base.xml

and

.. code-block:: console

  inputFiles/wellbore/CasedElasticWellbore_benchmark.xml

The corresponding integrated test is

.. code-block:: console

  inputFiles/wellbore/CasedElasticWellbore_smoke.xml

In this example, we would focus our attention on the ``Solvers``, ``Mesh`` and ``Constitutive`` tags.

-----------------------------------------------------------
Solid mechanics solver
-----------------------------------------------------------

As fluid flow is not considered, only the solid mechanics ``SolidMechanicsLagrangianSSLE`` solver is required for solving this linear elastic problem. In this solver, the three regions and three materials associated to casing, cement sheath and rock are respectively defined by ``targetRegions`` and ``solidMaterialNames``.  

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedElasticWellbore_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_SolidMechanicsSolver -->
  :end-before: <!-- SPHINX_SolidMechanicsSolverEnd -->


--------------------------------------------------------------------
Cased wellbore mesh
--------------------------------------------------------------------

The internal wellbore mesh generator ``InternalWellbore`` is employed to create the mesh of this wellbore problem. The radii of the casing cylinder, the cement sheath cylinder and the far-field boundary of the surrounding rock formation are defined by a vector ``radius``. In the tangent direction, ``theta`` angle is specified from 0 to 360 degree for a full geometry of the domain. Note that a half or a quarter of the domain can be defined by a ``theta`` angle from 0 to 180 or 90 degree, respectively. The trajectory of the well is defined by ``trajectory``, which is vertical in this case. The ``autoSpaceRadialElems`` parameters allow optimally increasing the element size from local zone around the wellbore to the far-field zone. In this example, the auto spacing option is only applied for the rock formation. The ``useCartesianOuterBoundary`` transforms the far-field boundary to a squared shape to enforce a Cartesian aligned outer boundary, which eases the loading of the boundary conditions. The ``cellBlockNames`` and ``elementTypes`` define the regions and related element types associated to casing, cement sheath and rock. 
 
.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedElasticWellbore_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_WellboreMesh -->
  :end-before: <!-- SPHINX_WellboreMeshEnd -->

.. figure:: mesh.png
   :align: center
   :figclass: align-center


--------------------------------------------------------------------   
Steel, cement, and rock constitutive laws
--------------------------------------------------------------------

Isotropic linear elastic constitutive behavior is considered for all the three materials. Note that the default density is useless for this case.

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_Material -->
  :end-before: <!-- SPHINX_MaterialEnd -->


--------------------------------------------------------------------	       
Boundary conditions
--------------------------------------------------------------------

Far-field boundary are subjected to roller constraints. The normal traction on the inner face of the casing is defined by ``Traction`` field specification. The nodeset generated by the internal wellbore generator for this face is named as ``rneg``. The traction type is ``normal`` to mimic a casing test pressure that is applied normal to the casing inner face . The negative sign of the scale is attributed to the negative sign convention for compressive stress in GEOS.

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_BoundaryConditions -->
  :end-before: <!-- SPHINX_BoundaryConditionsEnd -->

---------------------------------
Results and benchmark
---------------------------------

A good agreement between the GEOS results and analytical results is shown in the figure below:


.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/casedElasticWellbore/casedElasticWellboreFigure.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

This concludes the cased wellbore example.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
