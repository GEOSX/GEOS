.. _CompositionalMultiphaseFlow:

#######################################
Compositional Multiphase Flow Solver
#######################################

Introduction
=============

This flow solver is in charge of implementing the finite-volume discretization (mainly, accumulation and flux terms, boundary conditions) of the equations governing compositional multiphase flow in porous media.
The present solver can be combined with the :ref:`CompositionalMultiphaseWell` which handles the discrete multi-segment well model and provides source/sink terms for the fluid flow solver.

Below, we first review the set of :ref:`equations`, followed by a discussion of the
choice of :ref:`primary_variables` used in the global variable formulation.
Then we give an overview of the :ref:`discretization` and, finally, we provide a list of the solver :ref:`parameters` and an input :ref:`input_example`.

.. _theory:

Theory
=========================

.. _equations:

Governing Equations
-------------------

Mass Conservation Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mass conservation for component :math:`c` is expressed as:

.. math::
   \phi \frac{ \partial  }{\partial t} \bigg( \sum_\ell \rho_{\ell} \, y_{c \ell} \, S_{\ell} \bigg)
   + \nabla \cdot \bigg( \sum_\ell \rho_{\ell} \, y_{c \ell} \, \boldsymbol{u}_{\ell} \bigg)
   - \sum_\ell \rho_{\ell} \, y_{c \ell} \, q_{\ell} = 0,


where :math:`\phi` is the porosity of the medium,
:math:`S_{\ell}` is the saturation of phase :math:`\ell`, :math:`y_{c \ell}`
is the mass fraction of component :math:`c` in phase :math:`\ell`,
:math:`\rho_{\ell}` is the phase density, and :math:`t` is time. We note that the
formulation currently implemented in GEOS is isothermal.

Darcy's Law
~~~~~~~~~~~

Using the multiphase extension of Darcy's law, the phase velocity :math:`\boldsymbol{u}_{\ell}`
is written as a function of the phase potential gradient :math:`\nabla \Phi_{\ell}`:

.. math::
  \boldsymbol{u}_{\ell} := -\boldsymbol{k} \lambda_{\ell} \nabla \Phi_{\ell}
  = - \boldsymbol{k} \lambda_{\ell} \big( \nabla (p - P_{c,\ell}) - \rho_{\ell} g \nabla z \big).

In this equation, :math:`\boldsymbol{k}` is the rock permeability,
:math:`\lambda_{\ell} = k_{r \ell} / \mu_{\ell}` is the phase mobility,
defined as the phase relative permeability divided by the phase viscosity,
:math:`p` is the reference pressure, :math:`P_{c,\ell}` is the the capillary
pressure,  :math:`g` is the gravitational acceleration, and :math:`z` is depth.
The evaluation of the relative permeabilities, capillary pressures, and
viscosities is reviewed in the section about :doc:`/coreComponents/constitutive/docs/Constitutive`.

Combining the mass conservation equations with Darcy's law yields a set of :math:`n_c`
equations written as:

.. math::
   \phi \frac{ \partial  }{\partial t} \bigg( \sum_\ell \rho_{\ell} \, y_{c \ell} \, S_{\ell} \bigg)
   - \nabla \cdot \boldsymbol{k} \bigg( \sum_\ell \rho_{\ell} \, y_{c \ell} \, \lambda_{\ell} \nabla \Phi_{\ell}   \bigg)
   - \sum_\ell \rho_{\ell} \, y_{c \ell} \, q_{\ell} = 0.

Constraints and Thermodynamic Equilibrium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The volume constraint equation states that the pore space is always completely filled by
the phases. The constraint can be expressed as:

.. math::
  \sum_{\ell} S_{\ell} = 1.

The system is closed by the following thermodynamic equilibrium constraints:

.. math::
  f_{c \ell} - f_{c m} = 0.

where :math:`f_{c \ell}` is the fugacity of component :math:`c` in phase :math:`\ell`.
The flash calculations performed to enforce the thermodynamical equilibrium are reviewed
in the section about :doc:`/coreComponents/constitutive/docs/Constitutive`.

To summarize, the compositional multiphase flow solver assembles a set of :math:`n_c+1`
equations in each element, i.e., :math:`n_c` mass conservation equations and one volume constraint equation.
A separate module discussed in the :doc:`/coreComponents/constitutive/docs/Constitutive`
is responsible for the enforcement of the thermodynamic equilibrium at each nonlinear iteration.

==================== ===========================
Number of equations  Equation type
==================== ===========================
:math:`n_c`          Mass conservation equations
1                    Volume constraint
==================== ===========================

.. _primary_variables:

Primary Variables
------------------

The variable formulation implemented in GEOS is a global variable formulation based on
:math:`n_c+1` primary variables, namely, one pressure, :math:`p`, and
:math:`n_c` component densities, :math:`\rho_c`.
By default, we use molar component densities.
A flag discussed in the section :ref:`parameters` can be used to select mass component densities instead of molar component densities.

=========================== ===========================
Number of primary variables Variable type
=========================== ===========================
1                           Pressure
:math:`n_c`                 Component densities
=========================== ===========================

Assembling the residual equations and calling the
:doc:`/coreComponents/constitutive/docs/Constitutive` requires computing the molar component
fractions and saturations. This is done with the relationship:

.. math::
  z_c := \frac{\rho_c}{\rho_T},

where

.. math::
  \rho_T := \sum_c \rho_c.

These secondary variables are used as input to the flash calculations.
After the flash calculations, the saturations are computed as:

.. math::
  S_{\ell} := \nu_{\ell} \frac{ \rho_T }{ \rho_{\ell}},

where :math:`\nu_{\ell}` is the global mole fraction of phase :math:`\ell`
and :math:`\rho_{\ell}` is the molar density of phase :math:`\ell`.
These steps also involve computing the derivatives of the component
fractions and saturations with respect to the pressure and component densities.

.. _discretization:

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

The compositional multiphase solver uses a fully implicit (backward Euler) temporal discretization.

.. _solution_strategy:

Solution Strategy
-----------------

The nonlinear solution strategy is based on Newton's method.
At each Newton iteration, the solver assembles a residual vector, :math:`R`,
collecting the :math:`n_c` discrete mass conservation equations and the volume
constraint for all the control volumes.

.. _parameters:

Parameters
===========

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/CompositionalMultiphaseFVM.rst

.. _input_example:

Example
=========================

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
   :language: xml
   :start-after: <!-- START_SPHINX_INCLUDE_SOLVER_BLOCK -->
   :end-before: <!-- END_SPHINX_INCLUDE_SOLVER_BLOCK -->

We refer the reader to :ref:`TutorialDeadOilBottomLayersSPE10` for a complete tutorial illustrating the use of this solver.
