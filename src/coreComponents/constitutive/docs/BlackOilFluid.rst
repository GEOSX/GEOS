.. _BlackOilFluid:

############################################
Black-oil fluid model
############################################

Overview
=========================

In the black-oil model three pseudo-components, oil (o), gas (g) and water (w)
are considered. These are assumed to be partitioned across three fluid phases,
named liquid (l), vapor (v) and aqueous (a).

Phase behavior is characterized by the following quantities which are used to relate
properties of the fluids in the reservoir to their properties at surface conditions.

* :math:`B_o`: oil formation volume factor
* :math:`B_g`: gas formation volume factor
* :math:`R_s`: gas/oil ratio
* :math:`R_v`: oil/gas ratio

By  tables, that tabulate saturated and undersaturated oil and gas properties
as functions of pressure and solution ratios.

Dead oil
-------------
In **dead-oil** each component occupies only one phase. Thus, the following partition matrix determines the components distribution within the
three phases:

.. math::
    \begin{bmatrix}
    y_{gv} & y_{gl} & y_{ga}\\
    y_{ov} & y_{ol} & y_{oa}\\
    y_{wv} & y_{wl} & y_{wa}
    \end{bmatrix}
    = \begin{bmatrix}
    1 & 0 & 0 \\
    0  & 1 & 0 \\
    0 & 0 & 1
    \end{bmatrix}

and the phase densities are

.. math::
      \rho_{l} = & \, \frac{\rho_{o}^{STC}}{B_{o}} \\
      \rho_{v} = & \, \frac{\rho_{g}^{STC}}{B_{g}}.

Live oil
-------------
The live oil fluid model make no assumptions about the partitioning of the
hydrocarbon components and the following composition matrix can be used

.. math::
    \begin{bmatrix}
    y_{gv} & y_{gl} & y_{ga}\\
    \\
    y_{ov} & y_{ol} & y_{oa}\\
    \\
    y_{wv} & y_{wl} & y_{wa}
    \end{bmatrix}
    = \begin{bmatrix}
    \frac{\rho_{g}^{STC}}{\rho_{g}^{STC} + \rho_{o}^{STC} r_{s}} & \frac{\rho_{g}^{STC} R_{s}}{\rho_{o}^{STC} + \rho_{g}^{STC} R_{s}} & 0 \\
    \\
    \frac{\rho_{o}^{STC} r_{s}}{\rho_{g}^{STC} + \rho_{o}^{STC} r_{s}} & \frac{\rho_{o}^{STC}}{\rho_{o}^{STC} + \rho_{g}^{STC} R_{s}} & 0 \\
    \\
    0 & 0 & 1
    \end{bmatrix}

whereas the densities of the two hydrocarbon phases are

.. math::
  \rho_{l} = & \, \frac{\rho_{o}^{STC} + \rho_{g}^{STC} R_{s}}{B_{o}} \\
  \rho_{v} = & \, \frac{\rho_{g}^{STC} + \rho_{o}^{STC} R_{v}}{B_{g}}

See `Petrowiki`_ for more information.

Parameters
=========================

Both types are represented by ``<BlackOilFluid>`` node in the input.
Under the hood this is a wrapper around ``PVTPackage`` library, which is included as a submodule.
In order to use the model, GEOS must be built with ``-DENABLE_PVTPACKAGE=ON`` (default).

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/BlackOilFluid.rst

Supported phase names are:

===== ===========
Value Comment
===== ===========
oil   Oil phase
gas   Gas phase
water Water phase
===== ===========

Example
=========================

.. code-block:: xml

  <Constitutive>
    <BlackOilFluid name="fluid1"
                   fluidType="LiveOil"
                   phaseNames="{ oil, gas, water }"
                   surfaceDensities="{ 800.0, 0.9907, 1022.0 }"
                   componentMolarWeight="{ 114e-3, 16e-3, 18e-3 }"
                   tableFiles="{ pvto.txt, pvtg.txt, pvtw.txt }"/>
  </Constitutive>


.. _Petrowiki: https://petrowiki.spe.org/Phase_behavior_in_reservoir_simulation#Black-oil_PVT_models
