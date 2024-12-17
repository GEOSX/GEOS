Voigt Notation 
------------------------------------------------
In GEOS we express rank-two symmetric tensors using
`Voigt notation <https://en.wikipedia.org/wiki/Voigt_notation>`_.
Stress tensors are represented as an "unrolled" six-component
vector storing only the unique component values.  For strain tensors, note that *engineering
strains* are used such that the shear components of strain are multiplied by a
factor of two.
With this representation, the strain energy density is, conveniently, the inner product
of the stress and strain vectors.

Voigt representation of common rank-2 symmetric tensors are:

.. math::
   \bm{\sigma} = 
      \begin{bmatrix} \
         \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12}
      \end{bmatrix},
   \mathbf{S} = 
      \begin{bmatrix} 
         S_{11} \\ S_{22} \\ S_{33} \\ S_{23} \\ S_{13} \\ S_{12}
      \end{bmatrix},
   \bm{\epsilon} = 
      \begin{bmatrix} 
         \epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\\ 
         2 \epsilon_{23} \\ 2 \epsilon_{13} \\ 2 \epsilon_{12}
      \end{bmatrix}, 
   \mathbf{D} = 
      \begin{bmatrix}
         D_{11} \\ D_{22} \\ D_{33} \\ 2 D_{23} \\ 2 D_{13} \\ 2 D_{12}
      \end{bmatrix},
   \mathbf{E} = 
      \begin{bmatrix}
         E_{11} \\ E_{22} \\   E_{33} \\ 2 E_{23} \\ 2 E_{13} \\ 2 E_{12}
      \end{bmatrix},
   \mathbf{B} = 
      \begin{bmatrix}
         B_{11} \\ B_{22} \\ B_{33} \\ 2 B_{23} \\ 2 B_{13} \\ 2 B_{12}
   \end{bmatrix},
   \mathbf{C} = 
      \begin{bmatrix} 
        C_{11} \\ C_{22} \\ C_{33} \\ 2 C_{23} \\ 2 C_{13} \\ 2 C_{12}
      \end{bmatrix},

where
:math:`\bm{\sigma}` is the Cauchy stress,
:math:`\mathbf{S}` is the second Piola-Kirchhoff stress,
:math:`\bm{\epsilon}` is the small strain tensor,
:math:`\mathbf{D}` is the rate of deformation tensor,
:math:`\mathbf{E}` is the Lagrangian finite strain tensor,
:math:`\mathbf{B}` is the left Cauchy-Green tensor,
:math:`\mathbf{C}` is the right Cauchy-Green deformation tensor.

.. note::

   The factor of two in the shear components of strain (and strain-like) quantities is a frequent
   source of confusion, even for expert modelers.  It can be particularly challenging to use in nuanced 
   situations like stiffness tensor calculations or invariant decompositions.  If you plan to implement new models within 
   GEOS, please pay extra attention to this detail. We also provide many common operations in 
   centralized helper functions to avoid re-inventing the wheel.

