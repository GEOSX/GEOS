.. _BiotPorosity:


####################################################
Biot Porosity Model
####################################################


Overview
======================

According to the poroelasticity theory, the porosity (pore volume), :math:`\phi`, can be computed as

.. math::
   \phi = \phi_{ref} + \alpha ( \epsilon_v - \epsilon_{v,\,ref} ) + (p - p_{ref}) / N.

Here, :math:`\phi_{ref}` is the porosity at a reference state with pressure :math:`p_{ref}` and
volumetric strain :math:`\epsilon_{v,\,ref}`. Additionally, :math:`\alpha` is the Biot coefficient,
:math:`\epsilon_v` is the volumetric strain, :math:`p` is the fluid pressure
and :math:`N = \frac{K_s}{\alpha - \phi_{ref}}`, where :math:`{K_s}` is the grain bulk modulus.


Parameters
======================

The Biot Porosity Model can be called in the
``<Constitutive>`` block of the input XML file.
This porosity model must be assigned a unique name via the
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute in the ``<ElementRegions>``
block.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/BiotPorosity.rst

Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <BiotPorosity name="rockPorosity"
                    grainBulkModulus="1.0e27"
                    defaultReferencePorosity="0.3"/>
      ...
   </Constitutive>
