.. _CompressibleSinglePhaseFluid:

############################################
Compressible single phase fluid model
############################################

Overview
=========================

This model represents a compressible single-phase fluid with constant compressibility
and pressure-dependent viscosity.
These assumptions are valid for slightly compressible fluids, such as water, and some
types of oil with negligible amounts of dissolved gas.

Specifically, fluid density is computed as

.. math::
   \rho(p) = \rho_0 e^{c_\rho(p - p_0)}

where :math:`c_\rho` is compressibility, :math:`p_0` is reference pressure, :math:`\rho_0` is
density at reference pressure.
Similarly,

.. math::
   \mu(p) = \mu_0 e^{c_\mu(p - p_0)}

where :math:`c_\mu` is viscosibility (viscosity compressibility), :math:`\mu_0` is reference viscosity.

Either exponent may be approximated by linear (default) or quadratic terms of Taylor series expansion.
Currently there is no temperature dependence in the model, although it may be added in future.

Parameters
=========================

The model is represented by ``<CompressibleSinglePhaseFluid>`` node in the input.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/CompressibleSinglePhaseFluid.rst

Example
=========================

.. code-block:: xml

  <Constitutive>
    <CompressibleSinglePhaseFluid name="water"
                                  referencePressure="2.125e6"
                                  referenceDensity="1000"
                                  compressibility="1e-19"
                                  referenceViscosity="0.001"
                                  viscosibility="0.0"/>
  </Constitutive>
