.. _ThreePhaseRelativePermeability:

############################################
Three-phase relative permeability model
############################################


Overview
======================

For the simulation of three-phase flow in porous media, it is common to use a specific treatment
(i.e., different from the typical two-phase procedure) to evaluate the oil relative permeability.
Specifically, the three-phase oil relative permeability is obtained by interpolation of oil-water
and oil-gas experimental data measured independently in two-phase displacements.

Let :math:`k_{rw,wo}` and :math:`k_{ro,wo}` be the water-oil two-phase relative permeabilities for the
water phase and the oil phase, respectively. Let :math:`k_{rg,go}` and :math:`k_{ro,go}` be the oil-gas
two-phase relative permeabilities for the gas phase and the oil phase, respectively.
In the current implementation, the two-phase relative permeability data is computed analytically using the :doc:`/coreComponents/constitutive/docs/BrooksCoreyRelativePermeability`.

The water and gas three-phase relative permeabilities are simply given by two-phase data and
only depend on :math:`S_w` and :math:`S_g`, respectively. That is,

.. math::
    k_{rw,wog}(S_w) = k_{rw,wo}(S_w),

.. math::
    k_{rg,wog}(S_g) = k_{rg,go}(S_g).

The oil three-phase relative permeability
is obtained using a variant of the saturation-weighted interpolation procedure initially proposed
by `Baker <http://dx.doi.org/10.2118/17369-MS>`__. Specifically, we compute:

.. math::
    k_{ro,wog}(S_w,S_g) = \frac{ (S_w - S_{w,\textit{min}}) k_{ro,wo}(S_w) + S_g k_{rg,go}(S_g) }{ (S_w - S_{w,\textit{min}}) + S_g }.

This procedure provides a simple but effective formula avoiding
the problems associated with the other interpolation methods (negative values).

Another option can be triggered using `threePhaseInterpolator` to set interpolation model to be STONEII described by:

.. math::
    k_ro = k_rocw ((k_row/k_rocw + k_rw)(k_rog/k_rocw + k_rg) - k_rw - k_rg)

...

Parameters
======================

The relative permeability constitutive model is listed in the 
``<Constitutive>`` block of the input XML file.
The relative permeability model must be assigned a unique name via
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute of the ``<ElementRegion>``
node.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/BrooksCoreyBakerRelativePermeability.rst

Below are some comments on the model parameters.

* ``phaseNames`` - The number of phases should be 3. Supported phase names are:

===== ===========
Value Phase
===== ===========
oil   Oil phase
gas   Gas phase
water Water phase
===== ===========

* ``phaseMinVolFraction`` - The list of minimum volume fractions :math:`S_{\ell,min}` for each phase is specified in the same order as in ``phaseNames``. Below this volume fraction, the phase is assumed to be immobile.

* ``waterOilRelPermExponent`` - The list of exponents :math:`\lambda_{\ell,wo}` for the two-phase water-oil relative permeability data, with the water exponent first and the oil exponent next. These exponents are then used to compute :math:`k_{r \ell,wo}` in the :doc:`/coreComponents/constitutive/docs/BrooksCoreyRelativePermeability`.

* ``waterOilRelPermMaxValue`` - The list of maximum values :math:`k_{\textit{r} \ell,wo,\textit{max}}` for the two-phase water-oil relative permeability data, with the water max value first and the oil max value next. These exponents are then used to compute :math:`k_{r \ell,wo}` in the :doc:`/coreComponents/constitutive/docs/BrooksCoreyRelativePermeability`.

* ``gasOilRelPermExponent`` - The list of exponents :math:`\lambda_{\ell,go}` for the two-phase gas-oil relative permeability data, with the gas exponent first and the oil exponent next. These exponents are then used to compute :math:`k_{r \ell,go}` in the :doc:`/coreComponents/constitutive/docs/BrooksCoreyRelativePermeability`.

* ``gasOilRelPermMaxValue`` - The list of maximum values :math:`k_{\textit{r} \ell,go,\textit{max}}` for the two-phase gas-oil relative permeability data, with the gas max value first and the oil max value next. These exponents are then used to compute :math:`k_{r \ell,go}` in the :doc:`/coreComponents/constitutive/docs/BrooksCoreyRelativePermeability`.

Example
====================

.. code-block:: xml

   <Constitutive>
      ...
    <BrooksCoreyBakerRelativePermeability name="relperm"
                                          phaseNames="{oil, gas, water}"
                                          phaseMinVolumeFraction="{0.05, 0.05, 0.05}"
                                          waterOilRelPermExponent="{2.5, 1.5}"
                                          waterOilRelPermMaxValue="{0.8, 0.9}"
                                          gasOilRelPermExponent="{3, 3}"
                                          gasOilRelPermMaxValue="{0.4, 0.9}"/>
      ...
   </Constitutive>
