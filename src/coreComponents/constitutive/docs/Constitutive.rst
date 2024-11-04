.. _Constitutive:

Constitutive Models
============================================

Constitutive models describe relations between various physical quantities.
In a physics simulation they are used to model the response or state of material (solid, fluid, or a mixture) as a function of input variables.
There are many types of constitutive models available in GEOS.
These models are grouped together based on their input/output interface.

.. toctree::
   :maxdepth: 1

   solid/SolidModels
   FluidModels
   RelativePermeabilityModels
   CapillaryPressureModels
   PorosityModels
   PermeabilityModels
   PorousSolids
   TemperatureDependentSolidVolumetricHeatCapacity
   TemperatureDependentThermalConductivity


In an input XML file, constitutive models are listed in the ``<Constitutive>`` block.
Each parameterized model has its own XML tag, and each must be assigned a unique name via the ``name`` attribute.
Names are used to assign models to regions of the physical domain via the ``materialList`` attribute of the ``<CellElementRegion>`` node (see :ref:`XML_ElementRegions`).
In some cases, physics solvers must also be assigned specific constitutive models to use (see :ref:`Solvers`).

A typical ``<Constitutive>`` and ``<ElementRegions>`` block will look like:

.. code-block:: xml

  <Problem>

    <Constitutive>

      <!-- Define a compressible, single-phase fluid called "water"-->
      <CompressibleSinglePhaseFluid
        name="water"
        referencePressure="2.125e6"
        referenceDensity="1000"
        compressibility="1e-19"
        referenceViscosity="0.001"
        viscosibility="0.0"/>

    </Constitutive>

    <ElementRegions>

      <!--Add water to the material list for region 1-->
      <CellElementRegion
         name="region1"
         cellBlocks="{ * }"
         materialList="{ water }"/>

    </ElementRegions>

    ... remainder of problem definition here ...

  </Problem>
