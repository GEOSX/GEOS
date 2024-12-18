.. _ParallelPlatesPermeability:


####################################################
Parallel Plates Permeability Model
####################################################


Overview
======================

The parallel plates permeability model defines the relationship between the hydraulic
fracture aperture and its corresponding permeability following the classic lubrication model (`Witherspoon et al.`_ ) .
In this model, the two fracture walls are assumed to be smooth and parallel to each other and separated by a uniform aperture.

.. math::
   k =  \frac{a^3}{12}

where :math:`a` denotes the hydraulic fracture aperture.

Remark: :math:`k`, dimensionally, is not a permeability (as it is expressed in :math:`m^3`).



Parameters
======================

The Parallel Plates Permeability Model can be called in the
``<Constitutive>`` block of the input XML file.
This permeability model must be assigned a unique name via the
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute in the ``<ElementRegions>``
block.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/ParallelPlatesPermeability.rst


.. _Witherspoon et al.: https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/WR016i006p01016
