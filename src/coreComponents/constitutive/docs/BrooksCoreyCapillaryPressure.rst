.. _BrooksCoreyCapillaryPressure:

#####################################
Brooks-Corey capillary pressure model
#####################################


Overview
=====================

In GEOS, the oil-phase pressure is assumed to be the primary pressure.
The following paragraphs explain how the Brooks-Corey capillary pressure
model is used to compute the water-phase and gas-phase pressures as:

.. math::
    p_w = p_o - P_{c,w}(S_w),

and

.. math::
    p_g = p_o + P_{c,g}(S_g).

In the Brooks-Corey model, the water-phase capillary pressure
is computed as a function of the water-phase volume fraction with
the following expression:

.. math::
   P_{c,w}(S_w) = p_{e,w} S_{\textit{w,scaled}}^{-1/\lambda_w},

where the scaled water-phase volume fraction is computed as:

.. math::
   S_{\textit{w,scaled}} = \frac{S_w - S_{\textit{w,min}} }{1 - S_{\textit{w,min}} - S_{\textit{o,min}} - S_{\textit{g,min} }}.

The gas capillary pressure is computed analogously.

Parameters
====================

The capillary pressure constitutive model is listed in the 
``<Constitutive>`` block of the input XML file.
The capillary pressure model must be assigned a unique name via
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute of the ``<ElementRegions>``
node.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/BrooksCoreyCapillaryPressure.rst

Below are some comments on the model parameters:

* ``phaseNames`` - The number of phases can be either 2 or 3. The names entered for this attribute should match the phase names specified in the relative permeability block, either in :doc:`/coreComponents/constitutive/docs/BrooksCoreyRelativePermeability` or in :doc:`/coreComponents/constitutive/docs/ThreePhaseRelativePermeability`. The capillary pressure model assumes that oil is always present. Supported phase names are:

===== ===========
Value Phase
===== ===========
oil   Oil phase
gas   Gas phase
water Water phase
===== ===========

* ``phaseMinVolFraction`` - The list of minimum volume fractions :math:`S_{\ell,min}` for each phase is specified in the same order as in ``phaseNames``. Below this volume fraction, the phase is assumed to be immobile. The values entered for this attribute have to match those of the same attribute in the relative permeability block.

* ``phaseCapPressureExponentInv`` - The list of exponents :math:`\lambda_{\ell}` for each phase is specified in the same order as in ``phaseNames``. The parameter corresponding to the oil phase is currently not used.

* ``phaseEntryPressure`` - The list of entry pressures :math:`p_{e,\ell}` for each phase is specified in the same order as in ``phaseNames``. The parameter corresponding to the oil phase is currently not used.

* ``capPressureEpsilon`` - This parameter is used for both the water-phase and gas-phase capillary pressure. To avoid extremely large, or infinite, capillary pressure values, we set :math:`P_{c,w}(S_w) := P_{c,w}(\epsilon)` whenever :math:`S_w < \epsilon`. The gas-phase capillary pressure is treated analogously.

Example
======================

.. code-block:: xml

   <Constitutive>
      ...
      <BrooksCoreyCapillaryPressure name="capPressure"
                                    phaseNames="{oil, gas}"
                                    phaseMinVolumeFraction="{0.01, 0.015}"
                                    phaseCapPressureExponentInv="{0, 6}"
                                    phaseEntryPressure="{0, 1e8}"
                                    capPressureEpsilon="1e-8"/>
      ...
   </Constitutive>
