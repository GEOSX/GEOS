.. _ConstantPermeability:


####################################################
Constant Permeability Model
####################################################


Overview
======================

This model is used to define a diagonal permeability tensor that does not depend on any primary variable.

Parameters
======================


The following attributes are supported:

.. include:: /docs/sphinx/datastructure/ConstantPermeability.rst


Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <ConstantPermeability name="matrixPerm"
                            permeabilityComponents="{1.0e-12, 1.0e-12, 1.0e-12}"/>
      ...
   </Constitutive>
