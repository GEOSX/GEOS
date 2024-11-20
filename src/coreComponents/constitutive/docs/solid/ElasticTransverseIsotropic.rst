.. _ElasticTransverseIsotropic:

############################################
Model: Elastic Transverse-Isotropic
############################################

Overview
=========================

This model may be used for solid materials with a linear elastic, transverse-isotropic behavior.
This is most readily expressed in Voight notation as

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
      C_{11} & C_{12} &   C_{13} & 0 & 0 & 0 \\
      C_{12} & C_{11} &   C_{13} & 0 & 0 & 0 \\
      C_{13} & C_{13} &   C_{33} & 0 & 0 & 0 \\
      0 & 0 & 0 & C_{44} & 0 & 0 \\
      0 & 0 & 0 & 0 & C_{44} & 0 \\
      0 & 0 & 0 & 0 & 0 & (C_{11}-C_{12})/2
    \end{bmatrix}
    \begin{bmatrix}
      \epsilon_{11} \\
      \epsilon_{22} \\
      \epsilon_{33} \\
      2\epsilon_{23} \\
      2\epsilon_{13} \\
      2\epsilon_{12}
    \end{bmatrix}.

This system contains five independent constants.  These constants are calculated from the input parameters
indicated below.

Parameters
=========================

The following attributes are supported.  The "default" keyword in front of certain properties indicates that this
is the default value adopted for a region unless the user separately specifies a heterogeneous field via the
``FieldSpecification`` mechanism. 

.. include:: /docs/sphinx/datastructure/ElasticTransverseIsotropic.rst

Example
=========================

A typical ``Constititutive`` block will look like:

.. code-block:: xml

  <Constitutive>
    <ElasticTransverseIsotropic 
       name="shale"
       defaultDensity="2700"
       defaultPoissonRatioAxialTransverse="0.20"
       defaultPoissonRatioTransverse="0.30"
       defaultYoungModulusAxial="50.0e6"
       defaultYoungModulusTransverse="60.0e6"
       defaultShearModulusAxialTransverse="30.0e6" />
  </Constitutive>

