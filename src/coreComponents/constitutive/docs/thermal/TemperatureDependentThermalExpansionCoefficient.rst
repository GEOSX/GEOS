.. _TemperatureDependentThermalExpansionCoefficient:


##########################################################
Temperature-dependent Thermal Expansion Coefficient Model
##########################################################


Overview
======================

In this model, thermal expansion coefficient of porous medium is defined as a linear function of temperature: 

.. math::
   \alpha(T) = \alpha_{0} + \frac{ d \alpha }{ dT } (T - T_{0}) 

where 
      :math:`\alpha(T)` is the thermal expansion coefficient at temperature T; 
      :math:`\alpha_{0}` is the reference thermal expansion coefficient at the reference temperature; 
      :math:`T_{0}` is the reference temperature; 
      :math:`\frac{ d \alpha }{ dT }` is the gradient of the thermal expansion coefficient with respect to temperature, which equals to zero for the cases with constant thermal expansion coefficient.  


Parameters
======================

The temperature-dependent thermal expansion coefficient model is called in the
``<ElasticIsotropic>`` block of the input XML file.
This model must be assigned a unique name via the
``name`` attribute.
This name is used to attach the model to CellElementRegion of the physical
domain in the ``<ElementRegions>`` block.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/ElasticIsotropic.rst


Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <ElasticIsotropic
        name="rockSolid_linearTEC"
        defaultDensity="2700"
        defaultBulkModulus="9.12345679e9"
        defaultShearModulus="6.00813e9"
        defaultDrainedTEC="0.2e-5"
        dDrainedTEC_dT="0.0164e-5"
        referenceTemperature="0"/>
      ...
   </Constitutive>
