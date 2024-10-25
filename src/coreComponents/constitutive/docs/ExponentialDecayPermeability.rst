.. _ExponentialDecayPermeability:


####################################################
Exponential Decay Permeability Model
####################################################


Overview
======================

This stress-dependent permeability model assumes a simple exponential law for the fracture permeability as function of the effective normal stress acting on the fracture surface `(Gutierrez et al., 2000) <https://www.sciencedirect.com/science/article/pii/S0264817200000271>`__:

.. math::
   k =  k_{i} \text{exp} (- C {\sigma_n}')

where :math:`k_{i}` is the initial unstressed fracture permeability; :math:`C` is an empirical constant; :math:`{\sigma_n}'` is the effective normal stress.


Parameters
======================

The Exponential Decay Permeability model can be called in the
``<Constitutive>`` block of the input XML file.
This permeability model must be assigned a unique name via the
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``permeabilityModelName`` attribute in the ``<CompressibleSolidExponentialDecayPermeability>``
block.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/ExponentialDecayPermeability.rst


Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <ExponentialDecayPermeability 
        name="fracturePerm"
        empiricalConstant="0.27"
        initialPermeability="{1e-15, 1e-15, 1e-15}"/>
      ...
   </Constitutive>
