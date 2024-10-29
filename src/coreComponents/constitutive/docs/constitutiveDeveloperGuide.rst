.. _constitutiveModelsDoc:

################################################################################
Constitutive models
################################################################################

In GEOS, all constitutive models defining fluid and rock properties are implemented
in the namespace ``constitutive`` and derived from a common base class,
``ConstitutiveBase``. All objects are owned and handled by the ``ConstitutiveManager``.

Standalone models
======================================================
Standalone constitutive models implement constitutive laws such as:

  - mechanical material models (linear elasticity, plasticity, etc.),
  - PVT fluid behaviors,
  - relative permeability relationships,
  - porosity and permeability dependencies on state variables,
  - contact laws.

Storage, allocation, and update of properties
-----------------------------------------------------
Each constitutive model owns, as member variables, ``LvArray::Array`` containers
that hold the properties (or fields) and their derivatives with respect to the
other fields needed to update each property. Each property is stored as an array with the
first dimension representing the elementIndex and the second dimension storing the index of the
integration point. These dimensions are determined by the number of elements of the
subregion on which each constitutive model is registered, and by the chosen discretization
method. Vector and tensor fields have an additional dimension to identify
their components. Similarly, an additional dimension is
necessary for multiphase fluid models with properties defined for each component in each phase.
For example, a single-phase fluid model where density and viscosity are
functions of the fluid pressure has the following members:

.. literalinclude:: /coreComponents/constitutive/fluid/singlefluid/SingleFluidBase.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: //END_SPHINX_INCLUDE_00

Resizing all fields of the constitutive models happens during the initialization phase by
the ``ConstitutiveManger`` through a call to ``ConstitutiveManger::hangConstitutiveRelation``,
which sets the appropriate subRegion as the parent Group of each constitutive model object.
This function also resizes all fields based on the size of the subregion and the number of quadrature
points on it, by calling ``CONSTITUTIVE_MODEL::allocateConstitutiveData``. For the
single phase fluid example used before, this call is:

.. literalinclude:: /coreComponents/constitutive/fluid/singlefluid/SingleFluidBase.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: //END_SPHINX_INCLUDE_00

Any property or field stored on a constitutive model must be updated within a computational
kernel to ensure that `host` and `device` memory in GPUs are properly synced, and that any
updates are performed on `device`. Some properties are updated
within finite element kernels of specific physics (such as stress in a mechanics kernel). Consequently,
for each constitutive model class, a corresponding `nameOfTheModelUpdates`, which only contains
``LvArray::arrayView`` containers to the data, can be captured by value inside computational kernels.
For example, for the single phase fluid model `Updates` are:

.. literalinclude:: /coreComponents/constitutive/fluid/singlefluid/SingleFluidBase.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_01
   :end-before: //END_SPHINX_INCLUDE_01

Because `Updates` classes are responsible for updating the fields owned by the constitutive models,
they also implement all functions needed to perform property updates, such as:

.. literalinclude:: /coreComponents/constitutive/fluid/singlefluid/SingleFluidBase.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_02
   :end-before: //END_SPHINX_INCLUDE_02

Compound models
========================================================

Compound constitutive models are employed to mimic the behavior of a material that
requires a combination of constitutive models linked together. These compound models
do not hold any data. They serve only as an interface with the individual models that
they couple.

Coupled Solids
---------------------------------------------------------

``CoupledSolid`` models are employed to represent porous materials that require
both a mechanical behavior and constitutive laws that describe the
dependency of porosity and permeability on the primary unknowns.

The base class ``CoupledSolidBase`` implements some basic behaviors
and is used to access a generic ``CoupledSolid`` in a physics solver:

.. literalinclude:: /coreComponents/physicsSolvers/fluidFlow/SinglePhaseBase.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_COUPLEDSOLID
   :end-before: //END_SPHINX_INCLUDE_COUPLEDSOLID

Additionally, a `template class` defines a base ``CoupledSolid`` model
templated on the types of solid, porosity, and permeability models:

.. literalinclude:: /coreComponents/constitutive/solid/CoupledSolid.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: //END_SPHINX_INCLUDE_00

While physics solvers that need a porous material only interface with a compound model,
this one has access to the standalone models needed:

.. literalinclude:: /coreComponents/constitutive/solid/CoupledSolid.hpp
    :language: c++
    :start-after: //START_SPHINX_INCLUDE_01
    :end-before: //END_SPHINX_INCLUDE_01

There are two specializations of a ``CoupledSolid``:

- ``CompressibleSolid``: this model is used whenever there is no need to define a full mechanical model,
  but only simple correlations that compute material properties (like porosity or permeability).
  This model assumes that the solid model is of type `NullModel` and is only templated on
  the types of porosity and permeability models.

- ``PorousSolid``: this model is used to represent a full porous material where
  the porosity and permeability models need to be aware of the mechanical response of the material.
