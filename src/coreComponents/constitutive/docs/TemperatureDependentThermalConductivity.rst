.. _TemperatureDependentThermalConductivity:


##########################################################
Temperature-dependent Thermal Conductivity Model
##########################################################


Overview
======================

In this model, thermal conductivity of porous medium is defined as a linear function of temperature: 

.. math::
   k = k_{0} + \frac{ dk }{ dT } (T - T_{0}) 

where 
      :math:`k` is the thermal conductivity at temperature T, which is a vector; 
      :math:`k_{0}` is the reference thermal conductivity at the reference temperature; 
      :math:`T_{0}` is the reference temperature; 
      :math:`\frac{ dk }{ dT }` is the gradient of the thermal conductivity with respect to temperature, which equals to zero for the cases with constant thermal conductivity; note that this term is also in vector form, whose components could vary with directions.
  


Parameters
======================

The temperature-dependent thermal conductivity model is called in the
``<SinglePhaseThermalConductivity>`` block of the input XML file.
This model must be assigned a unique name via the
``name`` attribute.
This name is used to attach the model to CellElementRegion of the physical
domain in the ``<ElementRegions>`` block.

The following attributes are supported:

.. include:: /coreComponents/schema/docs/SinglePhaseThermalConductivity.rst


Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <SinglePhaseThermalConductivity
        name="thermalCond_nonLinear"
        defaultThermalConductivityComponents="{ 1.5, 1.5, 1.5 }"
        thermalConductivityGradientComponents="{ -12e-4, -12e-4, -12e-4 }"
        referenceTemperature="20"/>
      ...
   </Constitutive>
