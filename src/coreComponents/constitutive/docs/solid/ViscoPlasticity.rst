.. _ViscoPlasticModel:

############################################
Model: Viscoplasticity
############################################

The classical `Perzyna-type viscoplasticity models <https://www.sciencedirect.com/science/article/abs/pii/S0065215608700097>`__ are not suitable for rate-dependent viscoplastic models with nonsmooth multisurface, because of using unclearly defined nested viscoplastic loading surfaces. As an alternative, the `Duvaut-Lions viscoplastic theory <https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620261003>`__, precludes these difficulties by excluding the concept of nested viscoplastic loading surfaces. The viscoplastic constitutive equation that relates the stress :math:`\boldsymbol{\sigma}` and the viscoplastic strain rate :math:`\dot{\boldsymbol{\epsilon}^{vp}}` is given by:

.. math::
   \boldsymbol{\sigma} - \bar{\boldsymbol{\sigma}} = \frac{1}{t_*}\tensor{c}:\dot{\boldsymbol{\epsilon}^{vp}}

Here, :math:`\bar{\boldsymbol{\sigma}}` represents the inviscid stress, which is the rate-independent elasto-plastic stress part that can be solved by using elasto-plastic solvers (such as Drucker-Prager, CamClay, etc.). :math:`\tensor{c}` is the tangent stiffness tensor and :math:`t_*` is the relaxation time, which is measured in units of time. The viscoplastic strain rate :math:`\dot{\boldsymbol{\epsilon}}^{vp}` can be approximated using the following finite difference formula:

.. math::
   \dot{\boldsymbol{\epsilon}^{vp}} = \frac{1}{\Delta t}(\Delta \boldsymbol{\epsilon} - \Delta \boldsymbol{\epsilon}^{elas})

Here, :math:`\Delta t` is the time increment, :math:`\Delta \boldsymbol{\epsilon}` is the total strain increment, and :math:`\Delta \boldsymbol{\epsilon}^{elas}` is the elastic part of the strain increment. Note that the elastic strain increment is related to the stress increment through Hook's law.

.. math::
   \Delta \boldsymbol{\sigma} = \tensor{c}:\Delta \boldsymbol{\epsilon}^{elas}

With some arrangements, we can obtain the following formula to update the stress tensor of the Duvaut-Lions elasto-viscoplastic materials:

.. math::
   \boldsymbol{\sigma} = r_t \hat{\boldsymbol{\sigma}} + (1-r_t) \bar{\boldsymbol{\sigma}}

Here, the time ratio :math:`r_t` is calculated from the relaxation time :math:`t_*` and the time increment :math:`\Delta t` as:

.. math::
   r_t = \frac{1}{1+\Delta t/t_*}

Assuming elastic behavior, the trial stress tensor :math:`\hat{\boldsymbol{\sigma}}` is computed using the strain increment:

.. math::
   \hat{\boldsymbol{\sigma}}^{t+\Delta t} = \boldsymbol{\sigma}^t + \tensor{c}^{t+\Delta t}:\Delta \boldsymbol{\epsilon}^{t+\Delta t}

The tangent stiffness tensor is updated using the following equivalent approximation:

.. math::
   \tensor{c}^{t+\Delta t} = r_t \tensor{c}^e + (1-r_t) \tensor{c}^{t}

Here, :math:`\tensor{c}^e` is the elastic stiffness tensor.

The name of the viscoplastic solver is formed by adding the prefix `Visco` to the name of the elasto-plastic solver used to compute the inviscid stress :math:`\overline{\boldsymbol{\sigma}}`. For example, the solver `Visco Drucker-Prager` corresponds to cases where the inviscid stress is computed by the `Drucker-Prager` solver. It is interesting to note that equivalent viscoelastic solutions can also be obtained using the Duvault-Lions algorithm by updating the inviscid stress with an elastic solver.
