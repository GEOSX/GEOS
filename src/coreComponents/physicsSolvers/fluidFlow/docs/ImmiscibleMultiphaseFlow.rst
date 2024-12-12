.. _ImmiscibleMultiphaseFlow:

#######################################
Immiscible Multiphase Flow Solver
#######################################

Introduction
=============

This flow solver is used to implement the finite-volume discretization for the problem of modeling multiphase flow in porous media under the influence of viscous, gravity, and capillary forces while neglecting miscibility and taking into account rock and fluid compressibility.

In here, we go over the governing equations :ref:`immiscible_equations` that covers two different formulation options, followed by the :ref:`immiscible_discretization`, and we conclude by providing a list of the solver :ref:`immiscible_parameters` and an input :ref:`immiscible_input_example`.

.. _immiscible_theory:

Theory
=========================

.. _immiscible_equations:

Governing Equations
-------------------

Mass Conservation Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We consider a two-component system, say gas and water, flow in a compressible porous medium, in which both components can exist only in their corresponding phases of vapor and liquid. The gas and water components are denoted by the subscripts :math:`g` and
:math:`w`, respectively. Moreover, the liquid, which is the wetting phase, and the vapor, the non-wetting phase, are denoted by the subscripts :math:`\ell` and
:math:`v`, respectively. The mass conservation laws are expressed as:

.. math::
  \frac{\partial}{\partial t} (\phi\rho_v S_v) + \nabla \cdot (\rho_v \boldsymbol{u}_v) = 
  \rho_v q_v,

and

.. math::
  \frac{\partial}{\partial t} (\phi\rho_\ell S_\ell) + \nabla \cdot (\rho_\ell 
  \boldsymbol{u}_\ell) = \rho_\ell q_\ell,



where :math:`\phi(\mathbf{x})` is the porosity of the medium which is a function of pressure,
:math:`S_\ell(\mathbf{x},t)` is the saturation of the phase
:math:`\ell` and similarly for the phase :math:`v`, and :math:`t` is the time. The source/sink terms :math:`q_{\ell}` and :math:`q_{v}` are
positive for injection and negative for production. The phase
velocity, :math:`\boldsymbol{u}_\ell` and :math:`\boldsymbol{u}_v`, are defined using
the multiphase extension of Darcy's law (conservation of momentum) as

 .. math::
  \boldsymbol{u}_\ell := -k\lambda_\ell(\nabla p_\ell - \rho_\ell g \nabla z),

and

 .. math::
  \boldsymbol{u}_v := -k\lambda_v(\nabla p_v - \rho_v g \nabla z).

Here, :math:`k(\mathbf{x})` is the scalar absolute permeability of the medium, :math:`\lambda_\ell` is the phase mobility of the liquid phase defined as :math:`k_{r\ell}/\mu_\ell`, where :math:`k_{r\ell}(\mathbf{x},S_\ell)` is the phase relative permeability, :math:`\mu_\ell` is the phase viscosity, and :math:`\rho_{\ell}` is the phase density. 
These are also defined similarly for the vapor phase. In both cases we assume that the relative permeabilities are strictly increasing functions of their own saturation.
The gravitational acceleration is denoted by :math:`g`, and the
depth by :math:`z` (positive going downward).
The conservation of mass equations are constrained by the volume contraint equation:

.. math::
 S_{\ell} + S_v = 1,

Moreover, the capillary pressure constraint relates the two phase pressures with

.. math::
 P_{c}(S_{\ell}) = p_{v} - p_{\ell}.

We assume that capillary pressure is a strictly decreasing function of the wetting-phase saturation.

The evaluation of the relative permeabilities, capillary pressures, and
viscosities is reviewed in the section about :doc:`/coreComponents/constitutive/docs/Constitutive`.

We note that the formulation currently implemented in GEOS is isothermal. 

To summarize, the Immiscible multiphase flow solver assembles a set of :math:`n_p+1`
equations in each element, i.e., :math:`n_p` mass conservation equations and one volume constraint equation.

==================== ===========================
Number of equations  Equation type
==================== ===========================
:math:`n_p`          Mass conservation equations
1                    Volume constraint
==================== ===========================

.. _immiscible_primary_variables:

Primary Variables
------------------

There are two formulations implemented in GEOS for the Immiscible multiphsae solver and both formulations are based on
:math:`n_p+1` primary variables, namely, one pressure, :math:`p`, and
:math:`n_p` phase volume fractions, :math:`S_{p}`.


=========================== ===========================
Number of primary variables Variable type
=========================== ===========================
1                           Pressure
:math:`n_p - 1`             Phase volume fractions
=========================== ===========================

The main formulation is the standard formulation which solves the individual components mass conservation equations. Also, another formulation based on the total mass flux is implemented which is useful for multiple purposes such as hybrid upwinding techniques and sequential finite volume methods. This latter formulation is explained next. 

Flow and Transport Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To develop this formulation we use a flux approximation as required by the finite-volume numerical solution scheme.
Thus, we choose to construct this approximation in fractional flow form, and with this we will be able to show the coupling between the different physical processes. This formulation is obtained by decomposing the governing equations into a flow problem for both phases and a transport problem for one of the two phases. 
To obtain this decomposition,  we use a total-mass balance formulation by summing both components mass conversation equations and then using the mass constraint to result in the following elliptic PDE governing the temporal evolution of the pressure field:

.. math::
  \frac{\partial}{\partial t}(\phi \rho_t) + \nabla \cdot (\rho_{\ell} \boldsymbol{u}_{\ell} + \rho_{v} \boldsymbol{u}_{v}) = \rho_{\ell} q_{\ell} + \rho_{v} q_{v}, 

where 

.. math::
  \rho_t = \rho_{\ell}S_{\ell}+\rho_{v}S_{v}

and we defined a total mass flux as 

.. math::
 \boldsymbol{U}_T := \rho_{\ell} \boldsymbol{u}_{\ell} + \rho_{v} \boldsymbol{u}_{v}= -k (\rho_{\ell} \lambda_{\ell} + \rho_{v} \lambda_{v}) \nabla p + k ( \lambda_{\ell} \rho^2_{\ell} + \lambda_{v} \rho^2_{v}) g \nabla z + k \rho_{\ell}\lambda_{\ell} \nabla P_{c}. 

Next, the highly nonlinear parabolic transport equation is obtained by using  this total mass flux to formally eliminate the pressure variable from the individual components mass conservation equations, yielding

.. math::
 \frac{\partial}{\partial t}(\phi\rho_v S_v) + \nabla \cdot F_v
  = 
 \rho_v q_v,

and

.. math::

 \frac{\partial}{\partial t}(\phi\rho_\ell S_\ell) + \nabla \cdot F_\ell
  =  
 \rho_\ell q_\ell,

where the flow flux for each phase is defined as

.. math::
 {
 F_{\ell} :=
 \frac{\rho_\ell \lambda_\ell}{\rho_\ell \lambda_\ell+\rho_v 
 \lambda_v}\boldsymbol{U}_T}   + 
 k \frac{\rho_\ell \lambda_\ell\rho_v \lambda_v}{\rho_\ell \lambda_\ell+\rho_v 
 \lambda_v}(\rho_\ell - \rho_v) 
 g\nabla z
 +
 k \frac{\rho_\ell \lambda_\ell\rho_v \lambda_v}{\rho_\ell \lambda_\ell+\rho_v 
 \lambda_v} ( \nabla P_{c})

and

.. math::
 {
 F_{v} :=
 \frac{\rho_v \lambda_v}{\rho_\ell \lambda_\ell+\rho_v 
 \lambda_v}\boldsymbol{U}_T}   + 
 k \frac{\rho_\ell \lambda_\ell\rho_v \lambda_v}{\rho_\ell \lambda_\ell+\rho_v 
 \lambda_v}(\rho_v - \rho_\ell) 
 g\nabla z
 -
 k \frac{\rho_\ell \lambda_\ell\rho_v \lambda_v}{\rho_\ell \lambda_\ell+\rho_v 
 \lambda_v} ( \nabla P_{c})




.. _immiscible_discretization:

Discretization
--------------

Spatial Discretization
~~~~~~~~~~~~~~~~~~~~~~

The governing equations are discretized using standard cell-centered finite-volume
discretization.

In the approximation of the flux term at the interface between two control volumes,
the calculation of the pressure stencil is general and will ultimately support a
Multi-Point Flux Approximation (MPFA) approach. The current implementation of the
transmissibility calculation is reviewed in the section about
:doc:`/coreComponents/discretizationMethods/docs/NumericalMethodsManager`.

The approximation of the dynamic transport coefficients multiplying the discrete
potential difference (e.g., the phase mobilities) is performed with a first-order
phase-per-phase single-point upwinding based on the sign of the phase potential difference
at the interface.

Temporal Discretization
~~~~~~~~~~~~~~~~~~~~~~~

The immiscible multiphase solver uses a fully implicit (backward Euler) temporal discretization.

.. _immiscible_solution_strategy:

Solution Strategy
-----------------

The nonlinear solution strategy is based on Newton's method.
At each Newton iteration, the solver assembles a residual vector, :math:`R`,
collecting the :math:`n_p` discrete mass conservation equations and the volume
constraint for all the control volumes.

.. _immiscible_parameters:

Parameters
===========

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/ImmiscibleMultiphaseFlow.rst

.. _immiscible_input_example:

Example
=========================

.. literalinclude:: ../../../../../inputFiles/immiscibleMultiphaseFlow/immiscible_2phaseFlow_1d.xml
   :language: xml
   :start-after: <!-- START_SPHINX_INCLUDE_SOLVER_BLOCK -->
   :end-before: <!-- END_SPHINX_INCLUDE_SOLVER_BLOCK -->

