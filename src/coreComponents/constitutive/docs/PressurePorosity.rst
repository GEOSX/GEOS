.. _PressurePorosity:


####################################################
Pressure dependent porosity
####################################################


Overview
======================

This model assumes a simple exponential law for the porosity as function of pressure, i.e.

.. math::
   \phi = \phi_{ref} \, exp ( c \cdot ( p - p_{ref} ) )

where :math:`\phi_{ref}` is the reference porosity at  reference pressure, :math:`p_{ref}` , :math:`p` is the pressure and, :math:`c` is the compressibility.


Parameters
======================


The following attributes are supported:

.. include:: /docs/sphinx/datastructure/PressurePorosity.rst

Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <PressurePorosity name="rockPorosity"
                        referencePressure="1.0e27"
                        defaultReferencePorosity="0.3"
                        compressibility="1.0e-9"/>
      ...
   </Constitutive>
