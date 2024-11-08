.. _AdvancedExampleNonLinearThermalOedometric:


######################################################################
Non-linear thermal oedometric example with temperature dependent TEC
######################################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example uses the temperature-dependent Thermal Expansion Coefficient (TEC) constitutive model to solve a non-linear thermoelastic oedometric problem. Here, thermal convection and fluid flow are neglected by setting fluid pressure and fluid heat capacity to zero. A uniform temperature is applied with a linear increase from 0 to 90 °C over a period of time from 0 to :math:`10^5` (s).

The confining stress variation :math:`dP_{conf}` with respect to temperature change :math:`dT` is defined by the following tangent relationship:

.. math::
   dP_{conf} = \frac{E}{1-\nu} \alpha(T) dT

where :math:`E` and :math:`\nu` are the elastic Young's modulus and Poisson's ratio. :math:`\alpha(T)` is the temperature-dependent thermal expansion coefficient (TEC). 

Let us consider the particular case of a linear TEC: 

.. math::
   \alpha(T) = \alpha_0 + A (T - T_0) 

where :math:`\alpha_0` is the TEC at a reference temperature :math:`T_0` and :math:`A` is the constant gradient of TEC with respect to the temperature. The confining stress is obtained by integrating the equation above:

.. math::
   P_{conf} = \frac{E}{1-\nu} \frac{\alpha_0 + \alpha(T)}{2} (T - T_0)

This analytical result is used here to validate the numerical results obtained with GEOS.


**Input file**

This benchmark example uses no external input file and everything required is
contained within two GEOS xml files that are located at:

.. code-block:: console

  inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_base.xml

and

.. code-block:: console

  inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_benchmark.xml

The corresponding smoke test is

.. code-block:: console

  inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_smoke.xml

In this example, we explain the ``Constitutive`` and ``FieldSpecifications`` tags.

-----------------------------------------------------------
Constitutive
-----------------------------------------------------------

The thermoelastic properties are defined in the following xml block:  

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_base.xml
  :language: xml
  :start-after: <!-- SPHINX_ThermoElasticProperties -->
  :end-before: <!-- SPHINX_ThermoElasticPropertiesEnd -->

--------------------------------------------------------------------   
FieldSpecifications
--------------------------------------------------------------------

The uniform temperature applied to the system is defined in the following xml block:

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_base.xml
  :language: xml
  :start-after: <!-- SPHINX_TemperatureBC -->
  :end-before: <!-- SPHINX_TemperatureBCEnd -->

Where the temperature table is specified by

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_base.xml
  :language: xml
  :start-after: <!-- SPHINX_TemperatureTable -->
  :end-before: <!-- SPHINX_TemperatureTableEnd -->

---------------------------------
Results and benchmark
---------------------------------

A good agreement between the GEOS results and analytical results is shown in the figure below:


.. plot:: docs/sphinx/advancedExamples/validationStudies/thermoPoromechanics/nonLinearThermalOedometric/thermoElasticOedometricFigure.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

This concludes the example of pure thermal diffusion problem around a wellbore.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
