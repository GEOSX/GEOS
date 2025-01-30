.. _ContactMechanicsSolver:

#####################################
Contact Mechanics Solver
#####################################

Governing Equations
--------------------------

GEOS contact solvers solve the the balance of linear momentum within a fractured solid, accounting for the continuity of stress across surfaces (i.e., fractures), i.e.

.. math::

   \nabla \cdot \sigma = 0 \\\\
   [[\sigma]] \cdot \mathbf{n} = 0

Where:

* :math:`\sigma` is the stress tensor in the solid,
* :math:`\mathbf{n}` is the outward unit normal to the surface,
* :math:`[[\sigma]]` is the stress jump across the surface.

On each fracture surface, a no-interpenetration constraint is enforced. Additionally, tangential tractions can also be generated, which are modeled using a regularized Coulomb model to describe frictional sliding.

Solvers
--------------------------

There exist two broad classes of discretization methods that model fractures as lower dimensional entities 
(eg, 2D surfaces in a 3D domain): conforming grid methods and nonconforming (or embedded) methods. Both approaches have 
been implemented in GEOS in the following solvers:

.. toctree::
   :maxdepth: 1

   SolidMechanicsConformingFractures
   
   SolidMechanicsEmbeddedFractures
