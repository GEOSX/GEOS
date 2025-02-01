.. _TwoPhaseFluid:

############################################
Two-phase fluid model
############################################

Overview
=========================

This model represents a two-phase fluid with pressure-dependent density and viscosity.

For each phase, both density and viscosity are described as tabulated data, either in the form of ``TableFunction`` or text files.

In the case of text files, one file is expected per phase and should consist of three columns: pressure, density and viscosity.

Note that currently, there is no temperature dependence in the model.


Parameters
=========================

The model is represented by ``<TwoPhaseFluid>`` node in the input.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/TwoPhaseFluid.rst
              

Example using TableFunctions
============================

.. code-block:: xml

  <Constitutive>
    <TwoPhaseFluid
      name="fluid"
      phaseNames="{ oil, water }"
      densityTableNames="{ densityTableOil, densityTableWater }"
      viscosityTableNames="{ viscosityTableOil, viscosityTableWater }" />
  </Constitutive>

  <Functions>
    <TableFunction
      name="densityTableOil"
      coordinateFiles="{ pres_pvdo.txt }"
      voxelFile="dens_pvdo.txt"
      interpolation="linear" />

    <TableFunction
      name="viscosityTableOil"
      coordinateFiles="{ pres_pvdo.txt }"
      voxelFile="visc_pvdo.txt"
      interpolation="linear" />

    <TableFunction
      name="densityTableWater"
      coordinates="{ 2068000, 5516000, 30600000, 55160000 }"
      values="{ 980.683, 982.07, 992.233, 1002.265 }"
      interpolation="linear" />

    <TableFunction
      name="viscosityTableWater"
      coordinates="{ 0 }"
      values="{ 0.0003 }"
      interpolation="linear" />
  </Functions>  


Example using text files
=========================

.. code-block:: xml

  <Constitutive>
    <TwoPhaseFluid
      name="fluid"
      phaseNames="{ oil, water }"
      tableNames="{ oil.txt, water.txt }" />
  </Constitutive>


with, for example, ``water.txt`` being set as:

.. code-block:: text

  #  P(Pa) Dens(kg/m3) Visc(Pa.s)
   2068000     980.683     0.0003     
   5516000      982.07     0.0003
  30600000     992.233     0.0003
  55160000    1002.265     0.0003
