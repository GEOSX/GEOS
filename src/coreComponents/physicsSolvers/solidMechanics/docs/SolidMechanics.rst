.. _SolidMechanicsLagrangianFEM:

#####################################
Solid Mechanics Solver
#####################################

List of Symbols
===================

.. math::
   i,j,k &\equiv \text {indices over spatial dimensions} \notag \\
   a,b,c &\equiv \text {indices over nodes} \notag \\
   l &\equiv \text {indices over volumetric elements} \notag \\
   q,r,s &\equiv \text {indices over faces} \notag \\
   n &\equiv \text {indices over time} \notag \\
   kiter &\equiv \text {iteration count for non-linear solution scheme} \notag \\
   \Omega &\equiv \text {Volume of continuum body} \notag \\
   \Omega_{crack} &\equiv \text {Volume of open crack} \notag \\
   \Gamma &\equiv \text {External surface of } \Omega \notag \\
   \Gamma_t &\equiv \text {External surface where tractions are applied} \notag \\
   \Gamma_u &\equiv \text {External surface where kinematics are specified} \notag \\
   \Gamma_{crack} &\equiv \text {entire surface of crack} \notag \\
   \Gamma_{cohesive} &\equiv \text {surface of crack subject to cohesive tractions} \notag \\
   \eta_0 &\equiv \text {set of all nodes} \notag \\
   \eta_f & \equiv \text {set of all nodes on flow mesh}  \notag \\
   m &\equiv \text{mass} \notag \\
   \kappa_k &\equiv \text{all elements connected to element k} \notag \\
   \phi & \equiv \text { porosity} \notag \\
   p_f & \equiv \text { fluid pressure} \notag \\
   \mathbf{u} & \equiv \text { displacement} \notag \\
   \mathbf{q} & \equiv \text { volumetric flow rate} \notag \\
   \mathbf{T} & \equiv \text { Cauchy stress} \notag \\
   \rho & \equiv \text { density in the current configuration} \notag \\
   \mathbf{x}& \equiv \text { current position} \notag \\
   \mathbf{w}&\equiv \text { aperture, or gap vector} \notag


Introduction
============
The `SolidMechanics_LagrangianFEM` solver applies a Continuous Galerkin finite element method to solve the linear momentum balance equation.
The primary variable is the displacement field which is discretized at the nodes.

Theory
=========================

Governing Equations
--------------------------

The `SolidMechanics_LagrangianFEM` solves the equations of motion as given by

.. math::
   T_{ij,j} + \rho(b_{i}-\ddot{x}_{i}) = 0,

which is a 3-dimensional expression for the well known expression of Newtons Second Law (:math:`F = m a`).
These equations of motion are discretized using the Finite Element Method,
which leads to a discrete set of residual equations:

.. math::
   (R_{solid})_{ai}=\int\limits_{\Gamma_t} \Phi_a t_i   dA  - \int\limits_\Omega \Phi_{a,j} T_{ij}   dV +\int\limits_\Omega \Phi_a \rho(b_{i}-\Phi_b\ddot{x}_{ib})  dV = 0

Quasi-Static Time Integration
-----------------------------
The Quasi-Static time integration option solves the equation of motion after removing the inertial term, which is expressed by

.. math::
   T_{ij,j} + \rho b_{i} = 0,

which is essentially a way to express the equation for static equilibrium (:math:`\Sigma F=0`).
Thus, selection of the Quasi-Static option will yield a solution where the sum of all forces at a given node is equal to zero.
The resulting finite element discretized set of residual equations are expressed as

.. math::
   (R_{solid})_{ai}=\int\limits_{\Gamma_t} \Phi_a t_i   dA  - \int\limits_\Omega \Phi_{a,j} T_{ij}   dV + \int\limits_\Omega \Phi_a \rho b_{i}  dV = 0,

Taking the derivative of these residual equations wrt. the primary variable (displacement) yields

.. math::
    \pderiv{(R_{solid}^e)_{ai}}{u_{bj}} &=
            - \int\limits_{\Omega^e} \Phi_{a,k} \frac{\partial T_{ik}}{\partial u_{bj}}   dV,

And finally, the expression for the residual equation and derivative are used to express a non-linear system of equations

.. math::
   \left. \left(\pderiv{(R_{solid}^e)_{ai}}{u_{bj}} \right)\right|^{n+1}_{kiter}
   \left( \left. \left({u}_{bj} \right) \right|^{n+1}_{{kiter}+1} - \left. \left({u}_{bj} \right) \right|^{n+1}_{kiter} \right)
   = - (R_{solid})_{ai}|^{n+1}_{kiter} ,

which are solved via the solver package.

Implicit Dynamics Time Integration (Newmark Method)
---------------------------------------------------
For implicit dynamic time integration, we use an implementation of the classical Newmark method.
This update method can be posed in terms of a simple SDOF spring/dashpot/mass model.
In the following, :math:`M` represents the mass, :math:`C` represent the damping of the dashpot, :math:`K`
represents the spring stiffness, and :math:`F` represents some external load.

.. math::
   M a^{n+1} + C v^{n+1} + K u^{n+1} &= F_{n+1},  \\

and a series of update equations for the velocity and displacement at a point:

.. math::
   u^{n+1} &= u^n + v^{n+1/2} \Delta t,  \\
   u^{n+1} &= u^n + \left( v^{n} + \inv{2} \left[ (1-2\beta) a^n + 2\beta a^{n+1} \right] \Delta t \right) \Delta t, \\
   v^{n+1} &= v^n + \left[(1-\gamma) a^n + \gamma a^{n+1} \right] \Delta t.

As intermediate quantities we can form an estimate (predictor) for the end of step displacement and midstep velocity by
assuming zero end-of-step acceleration.

.. math::
   \tilde{u}^{n+1} &= u^n + \left( v^{n} + \inv{2}  (1-2\beta) a^n  \Delta t \right) \Delta t = u^n + \hat{\tilde{u}}\\
   \tilde{v}^{n+1} &= v^n + (1-\gamma) a^n  \Delta t =  v^n + \hat{\tilde{v}}

This gives the end of step displacement and velocity in terms of the predictor with a correction for the end
step acceleration.

.. math::
   u^{n+1} &= \tilde{u}^{n+1} + \beta a^{n+1} \Delta t^2 \\
   v^{n+1} &= \tilde{v}^{n+1} + \gamma a^{n+1} \Delta t

The acceleration and velocity may now be expressed in terms of displacement, and ultimately in terms
of the incremental displacement.

.. math::
   a^{n+1} &= \frac{1}{\beta \Delta t^2} \left(u^{n+1} - \tilde{u}^{n+1} \right)  = \frac{1}{\beta \Delta t^2} \left( \hat{u} - \hat{\tilde{u}} \right) \\
   v^{n+1} &= \tilde{v}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(u^{n+1} - \tilde{u}^{n+1} \right) = \tilde{v}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(\hat{u} - \hat{\tilde{u}} \right)

plugging these into equation of motion for the SDOF system gives:

.. math::
   M \left(\frac{1}{\beta \Delta t^2} \left(\hat{u} - \hat{\tilde{u}} \right)\right) + C \left( \tilde{v}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(\hat{u} - \hat{\tilde{u}} \right) \right) + K u^{n+1} &= F_{n+1}  \\

Finally, we assume Rayliegh damping for the dashpot.

.. math::
   C = a_{mass} M + a_{stiff} K

Of course we know that we intend to model a system of equations with many DOF.
Thus the representation for the mass, spring and dashpot can be replaced by our finite element discretized
equation of motion.
We may express the system in context of a nonlinear residual problem

.. math::
    (R_{solid}^e)_{ai} &=
        \int\limits_{\Gamma_t^e} \Phi_a t_i   dA  \\
        &- \int\limits_{\Omega^e} \Phi_{a,j} \left(T_{ij}^{n+1}+  a_{stiff} \left(\pderiv{T_{ij}^{n+1}}{\hat{u}_{bk}} \right)_{elastic} \left( \tilde{v}_{bk}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(\hat{u}_{bk} - \hat{\tilde{u}}_{bk} \right) \right) \right)  dV \notag \\
        &+\int\limits_{\Omega^e} \Phi_a \rho \left(b_{i}- \Phi_b  \left( a_{mass} \left( \tilde{v}_{bi}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(\hat{u}_{bi} - \hat{\tilde{u}}_{bi} \right) \right) + \frac{1}{\beta \Delta t^2}  \left( \hat{u}_{bi} - \hat{\tilde{u}}_{bi} \right) \right) \right)  dV ,\notag \\
    \pderiv{(R_{solid}^e)_{ai}}{\hat{u}_{bj}} &=
        - \int\limits_{\Omega^e} \Phi_{a,k} \left(\pderiv{T_{ik}^{n+1}}{\hat{u}_{bj}}+  a_{stiff} \frac{\gamma}{\beta \Delta t} \left(\pderiv{T_{ik}^{n+1}}{\hat{u}_{bj}} \right)_{elastic} \right)   dV \notag \\
        &- \left( \frac{\gamma a_{mass}}{\beta \Delta t} + \frac{1}{\beta \Delta t^2}  \right) \int\limits_{\Omega^e} \rho \Phi_a \Phi_c    \pderiv{ \hat{u}_{ci} }{\hat{u}_{bj}}dV .

Again, the expression for the residual equation and derivative are used to express a non-linear system of equations

.. math::
   \left. \left(\pderiv{(R_{solid}^e)_{ai}}{u_{bj}} \right)\right|^{n+1}_{kiter}
   \left( \left. \left({u}_{bj} \right) \right|^{n+1}_{{kiter}+1} - \left. \left({u}_{bj} \right) \right|^{n+1}_{kiter} \right)
   = - (R_{solid})_{ai}|^{n+1}_{kiter} ,

which are solved via the solver package. Note that the derivatives involving :math:`u` and :math:`\hat{u}` are interchangable,
as are differences between the non-linear iterations.

Explicit Dynamics Time Integration  (Special Implementation of Newmark Method with \gamma=0.5, \beta=0)
-------------------------------------------------------------------------------------------------------
For the Newmark Method, if \gamma=0.5, \beta=0, and the inertial term contains a diagonalized "mass matrix",
the update equations may be carried out without the solution of a system of equations.
In this case, the update equations simplify to a non-iterative update algorithm.

First the mid-step velocity and end-of-step displacements are calculated through the update equations

.. math::
   \tensor{v}^{n+1/2} &= \tensor{v}^{n} +  \tensor{a}^n \left( \frac{\Delta t}{2} \right), \text{ and} \\
   \tensor{u}^{n+1} &= \tensor{u}^n + \tensor{v}^{n+1/2} \Delta t.

Then the residual equation/s are calculated, and acceleration at the end-of-step is calculated via

.. math::
   \left( \tensor{M} + \frac{\Delta t}{2} \tensor{C} \right) \tensor{a}^{n+1} &=  \tensor{F}_{n+1} - \tensor{C} v^{n+1/2} - \tensor{K} u^{n+1} .

Note that the mass matrix must be diagonal, and damping term may not include the stiffness based damping
coefficient for this method, otherwise the above equation will require a system solve.
Finally, the end-of-step velocities are calculated from the end of step acceleration:

.. math::
   \tensor{v}^{n+1} &= \tensor{v}^{n+1/2} + \tensor{a}^{n+1} \left( \frac{\Delta t}{2} \right).

Note that the velocities may be stored at the midstep, resulting one less kinematic update.
This approach is typically referred to as the "Leapfrog" method.
However, in GEOS we do not offer this option since it can cause some confusion that results from the
storage of state at different points in time.


Parameters
=========================

In the preceding XML block, The `SolidMechanics_LagrangianFEM` is specified by the title of the subblock of the `Solvers` block.
The following attributes are supported in the input block for `SolidMechanics_LagrangianFEM`:

.. include:: /coreComponents/schema/docs/SolidMechanics_LagrangianFEM.rst

The following data are allocated and used by the solver:

.. include:: /coreComponents/schema/docs/SolidMechanics_LagrangianFEM_other.rst

Example
=========================

An example of a valid XML block is given here:

.. literalinclude:: ../../../../../inputFiles/solidMechanics/sedov_finiteStrain_smoke.xml
  :language: xml
  :start-after: <!-- SPHINX_SOLID_MECHANICS_SOLVER -->
  :end-before: <!-- SPHINX_SOLID_MECHANICS_SOLVER_END -->
