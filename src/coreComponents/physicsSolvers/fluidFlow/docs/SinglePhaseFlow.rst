.. _SinglePhaseFlow:

#####################################
Singlephase Flow Solver
#####################################

Introduction
=============

Here, we describe the single-phase flow solver.
The role of this solver is to implement the fully implicit finite-volume discretization (mainly, accumulation and source terms, boundary conditions) of the equations governing compressible single-phase flow in porous media.
This solver can be combined with the SinglePhaseWell class which handles the discrete multi-segment well model and provides source/sink terms for the fluid flow solver.


Theory
=========================

.. _singlephase_equations:

Governing Equations
-------------------

This is a cell-centered Finite Volume solver for compressible single-phase flow in porous media.
Fluid pressure as the primary solution variable.
Darcy's law is used to calculate fluid velocity from pressure gradient.
The solver currently only supports Dirichlet-type boundary conditions (BC) applied on cells or faces and Neumann no-flow type BC.

The following mass balance equation is solved in the domain:

.. math::
   \frac{\partial}{\partial t}(\phi\rho) + \boldsymbol{\nabla} \cdot (\rho\boldsymbol{u}) + q = 0,

where

.. math::
   \boldsymbol{u} = -\frac{1}{\mu}\boldsymbol{k}(\nabla p - \rho \boldsymbol{g})

and :math:`\phi` is porosity, :math:`\rho` is fluid density, :math:`\mu` is fluid viscosity,
:math:`\boldsymbol{k}` is the permeability tensor, :math:`\boldsymbol{g}` is the gravity vector,
and :math:`q` is the source function (currently not supported). The details on the computation of the density and the viscosity are given in :ref:`CompressibleSinglePhaseFluid`.

When the entire pore space is filled by a single phase, we can substitute the Darcy's law into the mass balance equation to obtain the single phase flow equation

.. math::
   \frac{\partial}{\partial t}(\phi\rho) - \boldsymbol{\nabla} \cdot \frac{\rho \boldsymbol{k}}{\mu} (\nabla p - \gamma \nabla z) + q = 0,

with :math:`\gamma \nabla z= \rho \boldsymbol{g}`.


.. _singlephase_discretization:

Discretization
--------------


Space Discretization
~~~~~~~~~~~~~~~~~~~~

Let :math:`\Omega \subset \mathbb{R}^n, \, n =1,2,3` be an open set defining the computational domain. We consider :math:`\Omega` meshed by element such that :math:`\Omega = \cup_{i}V_i` and integrate the single phase flow equation, described above, over each element :math:`V_i`:

.. math::
   \int_{V_i} \frac{\partial}{\partial t}(\phi\rho) dV - \int_{V_i} \boldsymbol{\nabla} \cdot \frac{\rho \boldsymbol{k}}{\mu} (\nabla p - \gamma \nabla z) dV + \int_{V_i} q dV  = 0.

Applying the divergence theorem to the second term leads to

.. math::
  \int_{V_i} \frac{\partial}{\partial t}(\phi\rho)_i - \oint_{S_i} \left(\frac{\rho \boldsymbol{k}}{\mu}(\nabla p -\gamma \nabla z)\right) \cdot \boldsymbol{n} dS + \int_{V_i} q dV  = 0.

where :math:`S_i` represents the surface area of the element :math:`V_i` and :math:`\boldsymbol{n}` is a outward unit vector normal to the surface.

For the flux term, the (static) transmissibility is currently computed with a Two-Point Flux Approximation (TPFA) as described in :ref:`FiniteVolume`.

The pressure-dependent mobility :math:`\lambda = \frac{\rho}{\mu}` at the interface is approximated using a first-order upwinding on the sign of the potential difference.

Time Discretization
~~~~~~~~~~~~~~~~~~~

Let :math:`t_0 < t_1 < \cdots < t_N=T` be a grid discretization of the time interval :math:`[t_0,T], \, t_0, T \in \mathbb{R}^+`. We use the backward Euler (fully implicit) method to integrate the single phase flow equation between two grid points :math:`t_n` and :math:`t_{n+1}, \, n< N` to obtain the residual equation:

.. math::
   \int_{V_i} \frac{(\phi\rho)_i^{n+1} - (\phi\rho)_i^n}{\Delta t} - \oint_{S_i} \left(\frac{\rho \boldsymbol{k}}{\mu}(\nabla p -\gamma \nabla z)\right)^{n+1} \cdot \boldsymbol{n} dS + \int_{V_i} q^{n+1} dV = 0

where :math:`\Delta t = t_{n+1}-t_n` is the time-step. The expression of this residual equation and its derivative are used to form a linear system, which is solved via the solver package.

Parameters
===================

The solver is enabled by adding a ``<SinglePhaseFVM>`` node in the Solvers section.
Like any solver, time stepping is driven by events, see :ref:`EventManager`.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/SinglePhaseFVM.rst

In particular:

* ``discretization`` must point to a Finite Volume flux approximation scheme defined in the Numerical Methods section of the input file (see :ref:`FiniteVolume`)
* ``fluidName`` must point to a single phase fluid model defined in the Constitutive section of the input file (see :ref:`Constitutive`)
* ``solidName`` must point to a solid mechanics model defined in the Constitutive section of the input file (see :ref:`Constitutive`)
* ``targetRegions`` is used to specify the regions on which the solver is applied

Primary solution field label is ``pressure``.
Initial conditions must be prescribed on this field in every region, and boundary conditions
must be prescribed on this field on cell or face sets of interest.


Example
=========================

.. literalinclude:: ../../../../../inputFiles/singlePhaseFlow/3D_10x10x10_compressible_base.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_SOLVERS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_SOLVERS_END -->

We refer the reader to :ref:`this page <TutorialSinglePhaseFlowWithInternalMesh>` for a complete tutorial illustrating the use of this solver.
