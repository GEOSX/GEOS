
.. _SolutionStrategy:

#####################################
Solution Strategy
#####################################

All  physics solvers share a common solution strategy for nonlinear time-dependent
problems. Here, we briefly describe the nonlinear solver and the timestepping
strategy employed.

Nonlinear Solver
=================================
At each time-step, the nonlinear system of discrete residual equations, i.e.

.. math::
    r(x) = 0

is solved by employing the Newton-Raphson method. Here, :math:`x` is the vector of
primary unknowns. Thus, each physics solver is responsible for assembling
the Jacobian matrix :math:`J` containing the analytical derivatives of the residual
vector :math:`r` with respect to the primary variables. Then, at each Newton iteration
:math:`\nu`, the following linear system is solved

.. math::
    J^{\nu} \delta x^{\nu+1} = -r^{\nu},

where, :math:`\delta x^{\nu+1}` is the Newton update. This linear system can be
solved with a variety of different linear solvers described in :doc:`/coreComponents/linearAlgebra/docs/LinearSolvers`.
The Newton update, :math:`\delta x^{\nu+1}` is then applied to the primary variables:

..  math::
  x^{\nu+1} = x^{\nu} + \delta x^{\nu+1}.

This procedure is repeated until convergence is achieved or until the maximum number of
iterations is reached.

Line Search
---------------------------
A line search method can be applied along with the Newton's method to facilitate Nonlinear
convergence. After the Newton update, if the residual norm has increased instead
of decreased, a line search algorithm is employed to correct the Newton update.


The user can choose between two different behaviors in case the line search fails
to provide a reduced residual norm:

1. accept the solution and move to the next Newton iteration;

2. reject the solution and request a timestep cut;

Timestepping Strategy
==================================

The actual timestep size employed is determined by a combination of several factors.
In particular, specific output events may have timestep requirements that force a
specific timestep to be used. However, physics solvers do have the possibility of
requesting a specific timestep size to the event manager based on their specific
requirements. In particular, in case of fast convergence indicated by a small number of
Newton iterations, i.e.

.. math::
     \text{numIterations} < \text{dtIncIterLimit} \cdot \text{newtonMaxIter},

the physics solver will require to double the timestep size. On the other hand,
if a large number of nonlinear iterations are necessary to
find the solution at timestep :math:`n`

.. math::
     \text{numIterations} > \text{dtCutIterLimit} \cdot \text{newtonMaxIter},

the physics solver will request the next timestep, :math:`n+1`, to be half the size of timestep :math:`n`.
Here,

Additionally, in case the nonlinear solver fails to converge with the timestep provided by the
event manager, the timestep size is cut, i.e.

.. math::
     \text{dt} = \text{timestepCutFactor} \cdot \text{dt},

and the nonlinear loop is repeated with the new timestep size.


Parameters
============================

All parameters defining the behavior of the nonlinear solver and determining the
timestep size requested by the physics solver are defined in the NonlinearSolverParameters
and are presented in the following table.

.. include:: /docs/sphinx/datastructure/NonlinearSolverParameters.rst
