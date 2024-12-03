.. _WillisRichardsPermeability:


####################################################
Willis-Richards Permeability Model
####################################################


Overview
======================

In the Willis-Richards permeability model, the stress-aperture relationship is derived based on Barton–Bandis constitutive model. In this model, fracture hydraulic aperture is assumed to be a function of effective normal stress acting on the fracture surface and shear displacement along the fracture surface `(Willis-Richards et al., 1996)  <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/96JB00882>`__. 

.. math::
   a = \frac{ a_{m} + U_{s} \text{tan} \left( \phi_{dil} \right) }{ 1 + 9 \frac{ \sigma_{n} }{ \sigma_{ref} }} 

Based on the assumption of parallel plates, the correlation between fracture hydraulic aperture and its corresponding permeability is defined as:

.. math::
   k =  \frac{a^2}{12} 

where 
      :math:`a` is the fracture hydraulic aperture; 
      :math:`a_{m}` is the fracture aperture at zero contact stress; 
      :math:`U_{s}` is the relative shear displacement; 
      :math:`\phi_{dil}` is the shear dilation angle; 
      :math:`\sigma_{n}` is the effective normal stress acting on the fracture surface; 
      :math:`\sigma_{ref}` is the effective normal stress that causes a 90% reduction in the fracture hydraulic aperture. 


Parameters
======================

The Willis-Richards permeability model is called in the
``<Constitutive>`` block of the input XML file.
This permeability model must be assigned a unique name via the
``name`` attribute.
This name is used to attach the model to regions of the physical
domain via a ``permeabilityModelName`` attribute in the ``<CompressibleSolidWillisRichardsPermeability>``
block.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/WillisRichardsPermeability.rst


Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <WillisRichardsPermeability
        name="fracturePerm"
        maxFracAperture="0.005"
        dilationCoefficient="0.01"
        refClosureStress="1.0e7"/>
      ...
   </Constitutive>
