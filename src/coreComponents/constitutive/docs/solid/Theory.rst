.. _DeformationTheory:

Deformation Theories
=======================

.. contents:: Table of Contents
    :depth: 3

Introduction
------------------

The solid mechanics solvers in GEOS work in a time-discrete setting, in which the system state
at time :math:`t^n` is fully known, and the goal of the solution procedure is to advance forward 
one timestep to :math:`t^{n+1} = t^n + \Delta t`.  
As part of this process, calls to a 
solid model must be made to compute the updated stress :math:`\bm{\sigma}^{n+1}` resulting from 
incremental deformation over the timestep.  
History-dependent models may also need to compute updates to one or more internal state 
variables :math:`Q^{n+1}`.

The exact nature of the incremental update will depend, however, on the kinematic
assumptions made. 
Appropriate measures of deformation and stress depend on assumptions of
`infinitesimal <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`_ or 
`finite <https://en.wikipedia.org/wiki/Finite_strain_theory>`_ 
strain, as well as other factors like rate-dependence and material anisotropy.

This section briefly reviews three main classes of solid models in GEOS, grouped by their kinematic assumptions. 
The presentation is deliberately brief, as much more extensive presentations can be 
found in almost any textbook on linear and nonlinear solid mechanics.


Small Strain Models
-------------------------------------------------

Let :math:`\bm{u}` denote the displacement field, and :math:`\nabla \bm{u}` its gradient. 
In small strain theory, ones assumes the displacement gradients :math:`\nabla \bm{u} \ll 1`.
In this case, it is sufficient to use the linearized strain tensor

.. math::

  \bm{\epsilon} = \frac{1}{2} \left( \nabla \bm{u} + \bm{u} \nabla \right )

as the deformation measure. Higher-order terms present in finite strain theories are neglected.
For inelastic problems, this strain is additively decomposed into elastic and inelastic components as

.. math::

  \bm{\epsilon} = \bm{\epsilon}^e + \bm{\epsilon}^{i}.

Inelastic strains can arise from a number of sources: plasticity, damage, etc.
Most constitutive models (including nonlinear elastic and inelastic models) can then be generically
expressed in rate form as

.. math::

  \dot{\bm{\sigma}} = \bm{c} : \dot{\bm{\epsilon}}^e

where :math:`\dot{\bm{\sigma}}` is the Cauchy stress rate and :math:`\bm{c}` is the tangent stiffness 
tensor.  Observe that the stress rate is driven by the elastic component :math:`\dot{\bm{\epsilon}}^e` of the strain rate.

In the time-discrete setting (as implemented in the code) the incremental constitutive update 
for stress is computed from a solid model update routine as

.. math::
   \bm{\sigma^{n+1}} = \bm{\sigma}(\Delta \bm{\epsilon}, \Delta t, Q^n),

where :math:`\Delta \bm{\epsilon} = \bm{\epsilon}^{n+1}-\bm{\epsilon}^n` is the incremental strain, 
:math:`\Delta t` is the timestep size (important for rate-dependent models), and
:math:`Q^n` is a collection of material state variables (which may include the previous stress and
strain).

For path and rate independent models, such as linear elasticity,
a simpler constitutive update may be formulated in terms of the total strain:

.. math::
   \bm{\sigma^{n+1}} = \bm{\sigma}(\bm{\epsilon^{n+1}}).

GEOS will use this latter form in specific, highly-optimized solvers when we know in advance that a
linear elastic model is being applied.  The more general interface is
the default, however, as it can accommodate a much wider range of constitutive behavior within a common
interface.

When implicit timestepping is used, the solid models must also provide the stiffness tensor,

.. math::
  \bm{c}^{n+1} = \frac{\partial \bm{\sigma}^{n+1}}{\partial \bm{\epsilon}^{n+1}},

in order to accurately linearize the governing equations.
In many works, this fourth-order tensor is referred to as the algorithmic or consistent tangent, in the
sense that it must be "consistent" with the discrete timestepping scheme being used
(`Simo and Hughes 1987 <https://doi.org/10.1016/0045-7825(85)90070-2>`_).  
For inelastic models, it depends not only on the intrinsic material stiffness, but also the incremental nature of the loading process.
The correct calculation of this stiffness can have a dramatic impact on the convergence rate of Newton-type
solvers used in the implicit solid mechanics solvers.

.. _DeformationTheory_Hypo:

Finite Deformation Models with Hypo-Materials
-------------------------------------------------

In the finite deformation regime, there are two broad classes of constitutive models frequently used:

- Hypo-elastic models (and inelastic extensions)
- Hyper-elastic models (and inelastic extensions)

Hypo-materials typically rely on a rate-form of the constitutive equations expressed in the spatial configuration.  
Let :math:`\bm{v}(\bm{x},t)` denote the spatial velocity field.  It can be decomposed into symmetric and anti-symmetric
components as

.. math::
   \bm{d} = \frac{1}{2} \left( \nabla \bm{v} + \bm{v} \nabla \right ) \qquad \text{and} \qquad 
   \bm{w} = \frac{1}{2} \left( \nabla \bm{v} - \bm{v} \nabla \right ),

where :math:`\bm{d}` is the deformation rate tensor and :math:`\bm{w}` is the spin tensor. 
A hypo-material model can be written in rate form as

.. math::
   \mathring{\bm{\tau}} = \bm{c} : \bm{d}^e

where :math:`\mathring{\bm{\tau}}` is an `objective rate <https://en.wikipedia.org/wiki/Objective_stress_rate>`_ of the Kirchoff stress 
tensor, :math:`\bm{c}` is the tangent stiffness tensor, 
and :math:`\bm{d}^e` is the elastic component of the deformation rate.
We see that the structure is similar to the rate form in the small strain regime, 
except the rate of Cauchy stress is replaced with an objective rate of Kirchoff stress, 
and the linearized strain rate is replaced with the deformation rate tensor.  
 
The key difference separating most hypo-models is the choice of the objective stress rate. 
In GEOS, we adopt the incrementally objective integration algorithm proposed by 
`Hughes and Winget (1980) <https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620151210>`__.
This method relies on the concept of an incrementally rotating frame of reference in order
to preserve objectivity of the stress rate. In particular, the stress update sequence is

.. math::

      \Delta{\tensor{R}} = ( \tensor{I} - \frac{1}{2} \Delta t {\tensor{w}} )^{-1} ( \tensor{I} + \frac{1}{2} \Delta t {\tensor{w}} )
      &\qquad \text{(compute incremental rotation)}, \\
      \tensor{\bar{\tau}}^{n} = \Delta{\tensor{R}} \tensor{\tau}^{n} \Delta{\tensor{R}}^T
      &\qquad \text{(rotate previous stress)}, \\
      \tensor{\tau}^{n+1} = \tensor{\bar{\tau}}^{n} + \Delta \tensor{\tau}
      &\qquad \text{(call constitutive model to update stress)}.

First, the previous timestep stress is rotated to reflect any rigid rotations occuring over the timestep.
If the model has tensor-valued state variables besides stress, these must also be rotated.
Then, a standard constitutive update routine can be called, typically driven by the incremental 
strain :math:`\Delta \bm{\epsilon} = \Delta t \bm{d}`.
In fact, an identical update routine as used for small strain models can be re-used at this point.

.. note::
   Hypo-models suffer from several well known
   deficiencies.  Most notably, the energy dissipation in a closed loading cycle of a hypo-elastic 
   material is not guaranteed to be zero, as one might desire from thermodynamic considerations.  

Finite Deformation Models with Hyper-Materials
-------------------------------------------------

Hyper-elastic models (and inelastic extensions) attempt to correct the thermodynamic deficiencies of their hypo-elastic cousins.
The constitutive update can be generically expressed at

.. math::
   \bm{S}^{n+1} = \bm{S}(\Delta \mathbf{F}, Q^n, \Delta t),

where :math:`\bm{S}` is the second Piola-Kirchoff stress and :math:`\Delta \mathbf{F}` is the incremental deformation gradient. 
Depending on the model, the deformation gradient can be converted to different deformation measures as needed.
Similarly, different stress tensors can be recovered through appropriate push-forward and pull-back operations.

In a hyperelastic material, the elastic response is 
expressed in terms of a stored strain-energy function that serves as the
potential for stress, e.g.

.. math::
   \mathbf{S} = \frac{\partial \psi (\tensor{C})}{ \partial \tensor{C} },

where :math:`\psi` is 
the stored energy potential, and :math:`\tensor{C}` is the right Cauchy-Green 
deformation tensor.  This potential guarantees that the energy dissipated or gained in a closed elastic cycle is zero.


