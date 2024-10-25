.. _KozenyCarmanPermeability:


####################################################
Kozeny-Carman Permeability Model
####################################################


Overview
======================

In the Kozeny-Carman model (see `ref`_), the permeability of a porous medium
is governed by several key parameters, including porosity, grain size, and grain shape:

.. math::
   k =  \frac{({s_{\epsilon}}{D_p})^2 {\phi}^3} {150({1-\phi})^2}

where :math:`s_{\epsilon}` is the sphericity of the particles, :math:`D_p` is the particle
diameter, :math:`\phi` is the porosity of the porous medium.


Parameters
======================

The Kozeny-Carman Permeability Model can be called in the
``<Constitutive>`` block of the input XML file.
This permeability model must be assigned a unique name via the
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute in the ``<ElementRegions>``
block.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/CarmanKozenyPermeability.rst



Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <CarmanKozenyPermeability name="matrixPerm"
                                particleDiameter="0.0002"
                                sphericity="1.0"/>
      ...
   </Constitutive>

.. _ref: https://en.wikipedia.org/wiki/Kozeny%E2%80%93Carman_equation
