.. _SolidMechanicsEmbeddedFractures:

#########################################
Solid mechanics embedded fractures solver
#########################################



Introduction
============


Discretization & soltuion strategy
==================================

The linear momentum balance equation is discretized using a low order finite element method. Moreover, to account for the influence of the fractures on the overall behavior, we utilize the enriched finite element method (EFEM)
with a piece-wise constant enrichment. This method employs an element-local enrichment of the FE space using the concept of assumedenhanced strain [1-6].


Example
=========================

An example of a valid XML block is given here:

.. literalinclude:: ../../../../../inputFiles/efemFractureMechanics/Sneddon_embeddedFrac_smoke.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_SOLVER -->
  :end-before: <!-- SPHINX_SNEDDON_SOLVER_END -->

Parameters
=========================

In the preceding XML block, The `SolidMechanicsEmbeddedFractures` is specified by the title of the subblock of the `Solvers` block. 
Note that the `SolidMechanicsEmbeddedFractures` always relies on the existance of a
The following attributes are supported in the input block for `SolidMechanicsEmbeddedFractures`:

.. include:: /docs/sphinx/datastructure/SolidMechanicsEmbeddedFractures.rst

The following data are allocated and used by the solver:

.. include:: /docs/sphinx/datastructure/SolidMechanicsEmbeddedFractures_other.rst

References
==========

1. Simo JC, Rifai MS. A class of mixed assumed strain methods and the method of incompatible modes. *Int J Numer Methods Eng.* 1990;29(8):1595-1638. Available at: http://arxiv.org/abs/https://onlinelibrary.wiley.com/doi/pdf/10.1002/nme.1620290802.

2. Foster CD, Borja RI, Regueiro RA. Embedded strong discontinuity finite elements for fractured geomaterials with variable friction. *Int J Numer Methods Eng.* 2007;72(5):549-581. Available at: http://arxiv.org/abs/https://onlinelibrary.wiley.com/doi/pdf/10.1002/nme.2020.

3. Wells G, Sluys L. Three-dimensional embedded discontinuity model for brittle fracture. *Int J Solids Struct.* 2001;38(5):897-913. Available at: https://doi.org/10.1016/S0020-7683(00)00029-9.

4. Oliver J, Huespe AE, Sánchez PJ. A comparative study on finite elements for capturing strong discontinuities: E-fem vs x-fem. *Comput Methods Appl Mech Eng.* 2006;195(37-40):4732-4752. Available at: https://doi.org/10.1002/nme.4814.

5. Borja RI. Assumed enhanced strain and the extended finite element methods: a unification of concepts. *Comput Methods Appl Mech Eng.* 2008;197(33):2789-2803. Available at: https://doi.org/10.1016/j.cma.2008.01.019.

6. Wu J-Y. Unified analysis of enriched finite elements for modeling cohesive cracks. *Comput Methods Appl Mech Eng.* 2011;200(45-46):3031-3050. Available at: https://doi.org/10.1016/j.cma.2011.05.008.
