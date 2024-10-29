.. _kgdViscosityDominated:


#######################################################
Viscosity dominated KGD hydraulic fracture
#######################################################

------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

The KGD problem addresses a single plane strain fracture growing in an infinite elastic domain. Basic assumptions and characteristic shape for this example is similar to those of another case (:ref:`kgdViscosityDominated`) except that the viscosity dominated regime is now considered. In this regime, more work is spent to move the fracturing fluid than to split the intact rock. In this test, slickwater with a constant viscosity of 1 cp is chosen as fracturing fluid, whose compressibility is neglected. To make the case identical to the viscosity dominated asymptotic solution, an ultra-low rock toughness :math:`K_{Ic}` is defined and fracture is assumed to be always propagating following fluid front. 
Asymptotic solutions of the fracture length :math:`\ell`, the net pressure :math:`p_0` and the fracture aperture :math:`w_0` at the injection point for the KGD fracture with a viscosity dominated regime are provided by `(Adachi and Detournay, 2002) <https://onlinelibrary.wiley.com/doi/abs/10.1002/nag.213?casa_token=Jp094LF1vgQAAAAA:cTFllaSy9ze1t6Cf2-PzUq2k51ZkM0dChlMRcsvkKzB1ILTxCboU4-0LrKv8Jao-Sx5t3O1hLiX38UEk>`__: 

.. math::
   \ell = 0.6152 (\frac{ E_p Q_0^3 }{ 12\mu })^{ 1/6 } t^{ 2/3 }

   w_0^2 = 2.1 (\frac{ 12\mu Q_0 }{ E_p })^{ 1/2 } \ell

   w_0 p_0 = 0.62 (12\mu Q_0 E_p)^{ 1/2 }

where the plane modulus :math:`E_p` is defined by

.. math:: E_p = \frac{ E }{ 1-\nu^2 }

and the term :math:`X` is given as:

.. math:: X = \frac{ 256 }{ 3 \pi^2 } \frac{ K_{Ic}^4 }{ \mu Q_0 {E_p}^3 }


**Input file**

The input xml files for this test case are located at:

.. code-block:: console

  inputFiles/hydraulicFracturing/kgdViscosityDominated_base.xml

and

.. code-block:: console

  inputFiles/hydraulicFracturing/kgdViscosityDominated_benchmark.xml

The corresponding integrated test with coarser mesh and smaller injection duration is also prepared:

.. code-block:: console

  inputFiles/hydraulicFracturing/kgdViscosityDominated_smoke.xml

Python scripts for post-processing and visualizing the simulation results are also prepared:

.. code-block:: console

  inputFiles/hydraulicFracturing/scripts/hydrofractureQueries.py

.. code-block:: console

  inputFiles/hydraulicFracturing/scripts/hydrofractureFigure.py


Fluid rheology and rock toughness are defined in the xml blocks below. Please note that setting an absolute zero value for the rock toughness could lead to instability issue. Therefore, a low value of :math:`K_{Ic}` is used in this example.

.. literalinclude:: ../../../../../../../inputFiles/hydraulicFracturing/kgdViscosityDominated_base.xml
  :language: xml
  :start-after: <!-- Sphinx_Solvers_SurfaceGenerator -->
  :end-before:  <!-- Sphinx_Solvers_SurfaceGenerator_End -->

.. literalinclude:: ../../../../../../../inputFiles/hydraulicFracturing/kgdViscosityDominated_base.xml
  :language: xml
  :start-after: <!-- Sphinx_Constitutive_CompressibleSinglePhaseFluid -->
  :end-before:  <!-- Sphinx_Constitutive_CompressibleSinglePhaseFluid_End -->


First, by running the query script 

.. code-block:: console

   python ./hydrofractureQueries.py kgdViscosityDominated

the HDF5 output is postprocessed and temporal evolution of fracture characterisctics (fluid pressure and fracture width at fluid inlet and fracure half length) are saved into a txt file ``model-results.txt``, which can be used for verification and visualization:

.. code-block:: console
		
  [['      time', '  pressure', '  aperture', '    length']]
           2 1.075e+06 0.0001176       1.5
           4 9.636e+05 0.0001645         2
           6 8.372e+05 0.0001917       2.5
           8  7.28e+05  0.000209         3
          10 6.512e+05  0.000222       3.5

Note: GEOS python tools ``geosx_xml_tools`` should be installed to run the query script (See `Python Tools Setup <https://geosx-geosx.readthedocs-hosted.com/projects/geosx-geospythonpackages/en/latest/>`_ for details). 

A good agreement between GEOS results and analytical solutions is shown in the comparison below, which is generated using the visualization script:

.. code-block:: console

   python ./kgdViscosityDominatedFigure.py


.. plot:: docs/sphinx/advancedExamples/validationStudies/hydraulicFracture/kgdViscosityDominated/kgdViscosityDominatedFigure.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

This concludes the viscosity dominated KGD example.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
