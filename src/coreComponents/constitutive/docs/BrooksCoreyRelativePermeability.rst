.. _BrooksCoreyRelativePermeability:

############################################
Brooks-Corey relative permeability model
############################################


Overview
======================

The following paragraphs explain how the Brooks-Corey
model is used to compute the phase relative permeabilities as a function
of volume fraction (i.e., saturation) with the expression:

.. math::
    k_{r\ell} = k_{\textit{r}\ell,\textit{max}} S_{\ell,\textit{scaled}}^{\lambda_{\ell}},

where the scaled volume fraction of phase :math:`\ell` is computed as:

.. math::
   S_{\ell,\textit{scaled}} = \frac{S_{\ell} - S_{\ell,\textit{min}} }{1 - \sum^{n_p}_{m=1} S_{\textit{m,min}} }.

The minimum phase volume fractions :math:`S_{\ell,\textit{min}}` are model parameters specified by the user.

Parameters
======================

The relative permeability constitutive model is listed in the
``<Constitutive>`` block of the input XML file.
The relative permeability model must be assigned a unique name via
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute of the ``<ElementRegions>``
node.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/BrooksCoreyRelativePermeability.rst

Below are some comments on the model parameters.

* ``phaseNames`` - The number of phases can be either two or three. Note that for three-phase flow, this model does not apply a special treatment to the intermediate phase relative permeability (no Stone or Baker interpolation). Supported phase names are:

===== ===========
Value Phase
===== ===========
oil   Oil phase
gas   Gas phase
water Water phase
===== ===========

* ``phaseMinVolFraction`` - The list of minimum volume fractions :math:`S_{\ell,min}` for each phase is specified in the same order as in ``phaseNames``. Below this volume fraction, the phase is assumed to be immobile.

* ``phaseRelPermExponent`` - The list of exponents :math:`\lambda_{\ell}` for each phase is specified in the same order as in ``phaseNames``.

* ``phaseMaxValue`` - The list of maximum values :math:`k_{\textit{r} \ell,\textit{max}}` for each phase is specified in the same order as in ``phaseNames``.


Examples
=======================

For a two-phase water-gas system (for instance in the CO2-brine fluid model), a typical relative permeability input looks like:

.. code-block:: xml

   <Constitutive>
      ...
      <BrooksCoreyRelativePermeability
        name="relPerm"
        phaseNames="{ water, gas }"
        phaseMinVolumeFraction="{ 0.02, 0.015 }"
        phaseRelPermExponent="{ 2, 2.5 }"
        phaseRelPermMaxValue="{ 0.8, 1.0 }"/>
      ...
   </Constitutive>

For a three-phase oil-water-gas system (for instance in the Black-Oil fluid model), a typical relative permeability input looks like:

.. code-block:: xml

   <Constitutive>
      ...
      <BrooksCoreyRelativePermeability
        name="relPerm"
        phaseNames="{ water, oil, gas }"
        phaseMinVolumeFraction="{ 0.02, 0.1, 0.015 }"
        phaseRelPermExponent="{ 2, 2, 2.5 }"
        phaseRelPermMaxValue="{ 0.8, 1.0, 1.0 }"/>
      ...
   </Constitutive>    
