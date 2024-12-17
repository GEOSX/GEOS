.. _Plasticity:

############################################
Plasticity Notation  
############################################

.. contents:: Table of Contents
    :depth: 3

Overview
--------

According to the `theory of plasticity  <https://en.wikipedia.org/wiki/Flow_plasticity_theory>`__ in the small strain regime, the total strain :math:`\boldsymbol{\epsilon}` can be additively split into elastic (:math:`\boldsymbol{\epsilon}^e`) and plastic (:math:`\boldsymbol{\epsilon}^p`) strains:

.. math::
   \boldsymbol{\epsilon} = \boldsymbol{\epsilon}^e + \boldsymbol{\epsilon}^p.

The plastic strain tensor is obtained from the flow rule: 

.. math::
   \dot{\boldsymbol{\epsilon}}^p=\dot{\lambda}\frac{\partial g}{\partial\boldsymbol{\sigma}},
   
in which :math:`\dot{\lambda} \geq 0` is the magnitude of plastic strain rate and :math:`g` is the plastic potential. The elastic strain is related to Cauchy stress tensor in rate form as:

.. math::
  \dot{\boldsymbol{\sigma}} = \tensor{c}^e : \dot{\boldsymbol{\epsilon}}^e,

where :math:`\tensor{c}^e` is the fourth order elastic stiffness tensor. The Cauchy stress tensor is related to the total strain as

.. math::
  \dot{\boldsymbol{\sigma}} = \boldsymbol{c}^{ep} : \dot{\boldsymbol{\epsilon}}, 
  
in which :math:`\tensor{c}^{ep}` is the fourth order elasto-plastic stiffness tensor.


Two-Invariant Models 
----------------------------------

Two-invariant plasticity models use the first invariant of the Cauchy stress tensor and the second invariant of the deviatoric stress tensor to describe the yield surface. 

Here we use the following stress invariants to define the yield surface:  the von Mises stress :math:`q = \sqrt{3J_2} = \sqrt{3/2} \|\boldsymbol{s}\|` and mean normal stress :math:`p = I_1/3`. Here, :math:`I_1` and :math:`J_2` are the first invariant of the stress tensor and second invariant of the deviatoric stress, defined as

.. math::
   I_1 = tr(\boldsymbol{\sigma})/3 \, , \quad J_2 = \frac{1}{2} \|\boldsymbol{s}\|^2 \, , \quad \boldsymbol{s}=\boldsymbol{\sigma}-p \boldsymbol{1} \, ,

in which :math:`\boldsymbol{1}` is the identity tensor. 

Similarly, we can define invariants of strain tensor, namely, volumetric strain :math:`\epsilon_v` and deviatoric strain :math:`\epsilon_s`.

.. math::
   \epsilon_v = tr(\boldsymbol{\epsilon}) \, , \quad   \epsilon_s = \sqrt{\frac{2}{3}} \| \boldsymbol{e}\|  \, , \, \quad \text{where} \, \quad \boldsymbol{e}=\boldsymbol{\epsilon}-\frac{1}{3} \epsilon_v \boldsymbol{1}.

Stress and strain tensors can then be recomposed from the invariants as:

.. math::
   \boldsymbol{\sigma} = p \, \boldsymbol{1} + \sqrt{\frac{2}{3}} q \, \hat{\boldsymbol{n}}

.. math::
   \boldsymbol{\epsilon} = \frac{1}{3} \epsilon_v \boldsymbol{1} + \sqrt{\frac{3}{2}}\epsilon_s \hat{\boldsymbol{n}}

in which :math:`\hat{\boldsymbol{n}} = \boldsymbol{e}/\|\boldsymbol{e}\|`.

The following two-invariant models are currently implemented in GEOS:

  - :ref:`DruckerPrager <DruckerPrager>`

  - :ref:`J2Plasticity <J2Plasticity>`

  - :ref:`ModifiedCamClay <ModifiedCamClay>`

  - :ref:`DelftEgg <DelftEgg>`

Three-Invariant Models
----------------------------------

Several three-invariant models are under active development, but are not yet available in develop.  If you are interested in helping to add additional material models, please submit a feature request.
