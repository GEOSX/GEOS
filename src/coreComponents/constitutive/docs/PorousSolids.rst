.. _PorousSolids:

############################################
Porous Solids
############################################

Overview
========================
Simulation of fluid flow in porous media and of poromechanics,
requires to define, along with fluid properties, the hydrodynamical properties of
the solid matrix. Thus, for porous media flow and and poromecanical simulation in GEOS,
two types of composite constitutive models can be defined to specify the characteristics
of a porous material: (1) a `CompressibleSolid` model, used for flow-only simulations and which
assumes that all poromechanical effects can be represented by the pressure dependency of the
porosity; (2) a `PorousSolid` model which, instead, allows to couple any solid model with
a `BiotPorosity` model and to include permeability's dependence on the mechanical response.


Both these composite  models require the names of the solid, porosity and permeability models
that, combined, define the porous material. The following sections outline how these models can be
defined in the Constitutive block of the xml input files and which type of submodels they
allow for.

CompressibleSolid
========================
This composite constitutive model requires to define a `NullModel` as solid model (since
no mechanical properties are used), a `PressurePorosity` model and any type of `Permeability` model.

To define this composite model the keyword `CompressibleSolid` has to be appended to the name
of the permeability model of choice, as shown in the following example for the `ConstantPermeability` model.


.. code-block:: xml

   <Constitutive>
     <CompressibleSolidConstantPermeability name="porousRock"
                                            solidModelName="nullSolid"
                                            porosityModelName="rockPorosity"
                                            permeabilityModelName="rockPermeability"/>

    <NullModel name="nullSolid"/>

    <PressurePorosity name="rockPorosity"
                      referencePressure="1.0e27"
                      defaultReferencePorosity="0.3"
                      compressibility="1.0e-9"/>

    <ConstantPermeability name="rockPermeability"
                          permeabilityComponents="{ 1.0e-4, 1.0e-4, 1.0e-4 }"/>

   </Constitutive>


PorousSolid
======================
To run poromechanical problems, the total stress is decomposed into an "effective stress" (driven by mechanical deformations) and a pore fluid
pressure component, following the `Biot theory of poroelasticity <https://doi.org/10.1016/B978-0-08-040615-2.50011-3>`__.
For single-phase flow, or multiphase problems with no capillarity, this decomposition reads

.. math::
   \sigma_{ij} = \sigma\prime_{ij}  - b p \delta_{ij}

where :math:`\sigma_{ij}` is the :math:`ij` component of the total stress tensor,
:math:`\sigma\prime_{ij}` is the :math:`ij` component of the effective (Cauchy) stress tensor,
:math:`b` is Biot's coefficient,
:math:`p` is fluid pressure,
and :math:`\delta` is the Kronecker delta.

The `PorousSolid` models simply append the keyword Porous in front of the solid model they contain,
e.g., PorousElasticIsotropic, PorousDruckerPrager, and so on. Additionally, they require to
define a `BiotPorosity` model and a `ConstantPermeability` model. For example, a Poroelastic material
with a certain permeability can be defined as

.. code-block:: xml

   <Constitutive>
     <PorousElasticIsotropic name="porousRock"
                             porosityModelName="rockPorosity"
                             solidModelName="rockSkeleton"
                             permeabilityModelName="rockPermeability"/>

     <ElasticIsotropic name="rockSkeleton"
                       defaultDensity="0"
                       defaultYoungModulus="1.0e4"
                       defaultPoissonRatio="0.2"/>

     <BiotPorosity name="rockPorosity"
                   grainBulkModulus="1.0e27"
                   defaultReferencePorosity="0.3"/>

     <ConstantPermeability name="rockPermeability"
                        permeabilityComponents="{ 1.0e-4, 1.0e-4, 1.0e-4 }"/>
   </Constitutive>

Note that any of the previously described solid models is used by the `PorousSolid` model
to compute the effective stress, leading to either poro-elastic, poro-plastic, or poro-damage
behavior depending on the specific model chosen.
