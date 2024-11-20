.. _TemperatureDependentSolidVolumetricHeatCapacity:


##########################################################
Temperature-dependent Solid Volumetric Heat Capacity Model
##########################################################


Overview
======================

In this model, solid volumetric specific heat capacity is assumed to be a linear function of temperature: 

.. math::
   C = C_{0} + \frac{ dC }{ dT } (T - T_{0}) 

where 
      :math:`C` is the solid volumetric heat capacity at temperature T; 
      :math:`C_{0}` is the reference solid volumetric heat capacity at the reference temperature; 
      :math:`T_{0}` is the reference temperature; 
      :math:`\frac{ dC }{ dT }` is the gradient of the volumetric heat capacity with respect to temperature, which equals to zero for the cases with constant solid volumetric heat capacity; 
  


Parameters
======================

The temperature-dependent solid volumetric heat capacity model is called in the
``<SolidInternalEnergy>`` block of the input XML file.
This model must be assigned a unique name via the
``name`` attribute.
This name is used to attach the model to regions of the physical
domain via a ``solidInternalEnergyModelName`` attribute in the ``<CompressibleSolidConstantPermeability>``
block.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/SolidInternalEnergy.rst



Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <SolidInternalEnergy
        name="rockInternalEnergy"
        referenceVolumetricHeatCapacity="4.56e6"
        dVolumetricHeatCapacity_dTemperature="1e6"
        referenceTemperature="0"
        referenceInternalEnergy="0"/>
      ...
   </Constitutive>
