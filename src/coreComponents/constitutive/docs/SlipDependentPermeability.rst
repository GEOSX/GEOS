.. _SlipDependentPermeability:


####################################################
Slip Dependent Permeability Model
####################################################


Overview
======================

The slip dependent permeability model defines the relationship between the relative shear displacement and fracture permeability. In this model, fractrues/faults are represented as slip interfaces.

.. math::
   k =  k_{i} \left[ (M_{mult} - 1) \text{tanh} \left( 3 \frac{ U_{s} }{ U_{s_{threshold}} } \right) + 1 \right] 

where :math:`k_{i}` is the initial fracture permeability; :math:`M_{mult}` is the maximum permeability multiplier; :math:`U_{s}` is the relative shear displacement; :math:`U_{s_{threshold}}` is the slip threshold.


Parameters
======================

The Slip Dependent Permeability model can be called in the
``<Constitutive>`` block of the input XML file.
This permeability model must be assigned a unique name via the
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``permeabilityModelName`` attribute in the ``<CompressibleSolidSlipDependentPermeability>``
block.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/SlipDependentPermeability.rst


Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <SlipDependentPermeability 
        name="fracturePerm"
        shearDispThreshold="0.005"
        maxPermMultiplier="1000.0"
        initialPermeability="{1e-15, 1e-15, 1e-15}"/>
      ...
   </Constitutive>
