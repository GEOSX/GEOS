.. _ElasticIsotropic:

############################################
Model: Elastic Isotropic
############################################

Overview
=========================

This model may be used for solid materials with a linear elastic isotropic behavior.
The relationship between stress and strain is given by `Hooke's Law <https://en.wikipedia.org/wiki/Hooke%27s_law>`__,
expressed as:

.. math::
   \sigma_{ij} = \lambda \epsilon_{kk} + 2 \mu \epsilon_{ij},

where :math:`\sigma_{ij}` is the :math:`ij` component of the Cauchy stress tensor,
:math:`\epsilon_{ij}` is the :math:`ij` component of the strain tensor,
:math:`\lambda` is the first Lamé elastic constant,
and :math:`\mu` is the elastic shear modulus.

Hooke's Law may also be expressed using `Voigt notation <https://en.wikipedia.org/wiki/Voigt_notation>`__ for stress and strain vectors as:

.. math::
   \tensor{\sigma} = \tensor{C} \cdot \tensor{\epsilon},

or,

.. math::
    \begin{bmatrix}
      \sigma_{11} \\
      \sigma_{22} \\
      \sigma_{33} \\
      \sigma_{23} \\
      \sigma_{13} \\
      \sigma_{12}
    \end{bmatrix}
    =
    \begin{bmatrix}
      2\mu+\lambda  &   \lambda     &   \lambda   & 0   & 0 & 0 \\
          \lambda     &  2\mu+\lambda   &   \lambda   & 0   & 0 & 0 \\
          \lambda     &    \lambda    & 2\mu+\lambda & 0  & 0 & 0 \\
          0         &       0     &       0 &\mu  & 0 & 0 \\
      0         &           0     & 0       & 0   & \mu & 0 \\
      0         &       0     & 0       & 0   & 0 & \mu
    \end{bmatrix}
    \begin{bmatrix}
      \epsilon_{11} \\
      \epsilon_{22} \\
      \epsilon_{33} \\
      2\epsilon_{23} \\
      2\epsilon_{13} \\
      2\epsilon_{12}
    \end{bmatrix}.

Variations
==========

For finite deformation solvers, the elastic isotropic model can be called within a hypo-elastic update routine. 
See :ref:`DeformationTheory_Hypo`

Parameters
=========================

The following attributes are supported.  Note that any two elastic constants can be provided, and the other
two will be internally calculated.  The "default" keyword in front of certain properties indicates that this
is the default value adopted for a region unless the user separately specifies a heterogeneous field via the
``FieldSpecification`` mechanism. 

.. include:: /docs/sphinx/datastructure/ElasticIsotropic.rst

Example
=========================

A typical ``Constititutive`` block will look like:

.. code-block:: xml

  <Constitutive>
    <ElasticIsotropic 
       name="shale"
       defaultDensity="2700"
       defaultBulkModulus="60.0e6"
       defaultShearModulus="30.0e6" />
  </Constitutive>
