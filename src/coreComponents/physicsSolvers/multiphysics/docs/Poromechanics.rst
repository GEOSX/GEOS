.. _PoroelasticSolver:

##########################
Poromechanics Solver
##########################

Introduction
===========================================

This section describes the use of the poroelasticity models implemented in GEOS.

Theory
=========================

Governing Equations
--------------------------

In our model, the geomechanics (elasticity) equation is expressed in terms of the total stress :math:`\mathbf{\sigma}`:

.. math::
    \nabla \mathbf{\sigma} + \rho_b \mathbf{g} = 0

where it relates to effective stress :math:`\mathbf{\sigma\prime}` and pore pressure :math:`p` through Biot's coefficient :math:`b`:

.. math::
    \mathbf{\sigma} = \mathbf{\sigma\prime} - b p\mathbf{I}

The fluid mass conservation equation is expressed in terms of pore pressure and volumetric (mean) total stress:

.. math::
    \left( \frac{1}{M} + \frac{b^2}{K_{dr}} \right) \frac{\partial p}{\partial t} + \frac{b}{K_{dr}} \frac{\partial \sigma_v}{\partial t} + \nabla \cdot \mathbf{v}_f = f

where :math:`M` is the Biot's modulus and :math:`K_{dr}` is the drained bulk modulus.

Unlike the conventional reservoir model that uses Lagranges porosity, in the coupled geomechanics and flow model, Eulers porosity :math:`\phi` is adopted so the porosity variation is derived as:

.. math::
    \partial \phi = \left( \frac{b-\phi}{K_s}\right) \partial p + \left( b-\phi \right) \partial \epsilon_v

where :math:`K_{s}` is the bulk modulus of the solid grain and :math:`\epsilon_v` is the volumetric strain.

Parameters
===========================================

The poroelasticity model is implemented as a main solver listed in
``<Solvers>`` block of the input XML file that calls both SolidMechanicsLagrangianSSLE and SinglePhaseFlow solvers.
In the main solver, it requires the specification of solidSolverName, flowSolverName, and couplingTypeOption.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/SinglePhasePoromechanics.rst

* ``couplingTypeOption``: defines the coupling scheme.

The solid constitutive model used here is PoroLinearElasticIsotropic, which derives from ElasticIsotropic and includes an additional parameter: Biot's coefficient. The fluid constitutive model is the same as SinglePhaseFlow solver. For the parameter setup of each individual solver, please refer to the guideline of the specific solver.

An example of a valid XML block for the constitutive model is given here:

.. literalinclude:: ../../../../../inputFiles/poromechanics/PoroElastic_Terzaghi_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_POROELASTIC_CONSTITUTIVE -->
  :end-before: <!-- SPHINX_POROELASTIC_CONSTITUTIVE_END -->

Example
===========================================

.. literalinclude:: ../../../../../inputFiles/poromechanics/PoroElastic_Terzaghi_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_POROELASTIC_SOLVER -->
  :end-before: <!-- SPHINX_POROELASTIC_SOLVER_END -->
