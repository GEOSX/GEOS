.. _AdvancedExampleViscoModifiedCamClay:


#############################################################################
Visco Modified CamClay model: Triaxial Driver versus Semi-Analytical Solution 
#############################################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example uses the Triaxial Driver to simulate a visco-elasto-plastic oedometric compression test of a Visco Modified CamClay solid. Oedometric condition with zero lateral strain together with loading/unloading axial strain periods are imposed. Semi-analytical results for the mean and shear stress variations :math:`\Delta p` and :math:`\Delta q` can be established, considering the Perzyna approach, for the imposed oedometric boundary conditions as `(Runesson et al. 1999) <https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1484(199901)4:1%3C75::AID-CFM60%3E3.0.CO;2-4>`__ :

.. math::
   \Delta p = K(\Delta\varepsilon_{V} - \Delta\lambda \frac{\partial G}{\partial p})

.. math::
   \Delta q = 2\mu(\Delta\varepsilon_{V} - \Delta\lambda \frac{3}{2}\frac{\partial G}{\partial q})

where :math:`K` and :math:`\mu` are elastic bulk and shear moduli, :math:`G` is the plastic potential and :math:`\Delta\lambda` is the visco-plastic multiplier that can be approximated by:

.. math::
   \Delta\lambda = \frac{\Delta t}{t_*} \frac{F}{3\mu\frac{\partial F}{\partial q}\frac{\partial G}{\partial q} + K\frac{\partial F}{\partial p}\frac{\partial G}{\partial p} + h}

in which :math:`\Delta t` is the time increment, :math:`t_*` is the relaxation time, :math:`F` is the stress function defining the visco-plastic yield surface and :math:`h` is the hardening rate defined by:

.. math::
   h = -\frac{\partial F}{\partial \lambda}


These solutions are implemented in a Python script associated to this example for verifying GEOS results.


**Input files**

This validation example uses two GEOS xml files that are located at:

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_base.xml

and

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_ViscoModifiedCamClay.xml

It also uses a set of table files located at:

.. code-block:: console

  inputFiles/triaxialDriver/tables/


A Python script for the semi-analytical solutions presented above as well as for post-processing the GEOS results is provided at:

.. code-block:: console

  src/docs/sphinx/advancedExamples/validationStudies/viscoplasticity/ViscoModifiedCamClay/TriaxialDriver_vs_SemiAnalytic_ViscoModifiedCamClay.py

For this example, we focus on the ``Task`` and the ``Constitutive`` tags.

------------------------------------------------------------------
Task
------------------------------------------------------------------

The imposed axial strain loading/unloading periods, the lateral zero strain, and the initial stress are defined in the ``Task`` block as:

.. literalinclude:: ../../../../../../../inputFiles/triaxialDriver/triaxialDriver_ViscoModifiedCamClay.xml
    :language: xml
    :start-after: <!-- SPHINX_TASK -->
    :end-before: <!-- SPHINX_TASK_END -->

------------------------------
Constitutive laws
------------------------------

The elasto-visco-plastic parameters are defined as:


.. literalinclude:: ../../../../../../../inputFiles/triaxialDriver/triaxialDriver_base.xml
    :language: xml
    :start-after: <!-- SPHINX_MATERIAL_VISCO_MODIFIED_CAMCLAY -->
    :end-before: <!-- SPHINX_MATERIAL_VISCO_MODIFIED_CAMCLAY_END -->


All constitutive parameters such as density, viscosity, and the bulk and shear moduli are specified in the International System of Units.

--------------------------------------------------------------
A comparison between GEOS results and semi-analytical results
--------------------------------------------------------------

The simulation results are saved in a text file, named ``ViscoModifiedCamClayResults.txt``. A comparison between the results given by the TriaxialDriver solver in GEOS and the semi-analytical results presented above is shown below. The discrepancy between these results may due to the difference between the Duvaut-Lions approach and the Perzyna approach for time dependant behavior when applying for the Modified CamClay model as discussed by `Runesson et al. (1999) <https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1484(199901)4:1%3C75::AID-CFM60%3E3.0.CO;2-4>`__. 


.. plot:: docs/sphinx/advancedExamples/validationStudies/viscoplasticity/ViscoModifiedCamClay/TriaxialDriver_vs_SemiAnalytic_ViscoModifiedCamClay.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
