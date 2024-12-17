.. _ExampleThermalLeakyWell:


#########################################################################
Non-isothermal CO2 Plume Evolution and Leakage Through an Abandoned Well
#########################################################################


**Context**

This validation case is a more complex version of the benchmark problem presented in :ref:`ExampleIsothermalLeakyWell`.
While the latter is based on simple isothermal and immiscible fluid properties, the present validation case
relies on a more realistic fluid behavior accounting for thermal effects and mass exchange between phases.
This non-isothermal benchmark test has been used in
`(Class et al., 2009) <https://link.springer.com/article/10.1007/s10596-009-9146-x>`__
to compare different implementations of CO2-brine fluid properties in the
context of CO2 injection and storage in saline aquifers.

Our goal is to review the sections of the XML file that are used to parameterize the CO2-brine fluid behavior,
and to demonstrate that GEOS produces similar results as those presented in
`(Class et al., 2009) <https://link.springer.com/article/10.1007/s10596-009-9146-x>`__.

**Input file**

This benchmark test is based on the XML file located below:

.. code-block:: console

  inputFiles/compositionalMultiphaseFlow/benchmarks/thermalLeakyWell/thermalLeakyWell_benchmark.xml

------------------------------------------------------------------------
Problem description
------------------------------------------------------------------------

Some of the text below is adapted from
`(Ebigbo, Class, Helmig, 2007) <https://link.springer.com/article/10.1007%2Fs10596-006-9033-7>`__.

The benchmark scenario remains the same as in :ref:`ExampleIsothermalLeakyWell`.
CO2 is injected into an aquifer, spreads within the aquifer, and, upon reaching a leaky well,
rises up to a shallower aquifer.
The model domain still has the dimensions: 1000 x 1000 x 160 m, but it is now assumed to be
shallower, between 640 m and 800 m of depth.

The figure below shows the pressure and temperature in the formation at the mentioned
depths (assuming a geothermal gradient of 0.03 K/m).
The conditions in the aquifer at the considered depths range from supercritical to liquid
to gaseous.
The figure also shows the CO2 density at the conditions of the formation.
There is a large change in density at a certain  depth.
This depth corresponds to the point where the line depicting the formation conditions
crosses the CO2 saturation vapor curve, that is, the boundary between liquid and gaseous CO2.
Other fluid properties such as viscosity also change abruptly at that depth.    

.. _thermalLeakyWell_aquiferConditions:
.. figure:: aquiferConditions.png
   :align: center
   :width: 500
   :figclass: align-center

   Aquifer conditions (image taken from `(Ebigbo, Class, Helmig, 2007) <https://link.springer.com/article/10.1007%2Fs10596-006-9033-7>`__).

Therefore, as explained later, we use a more sophisticated fluid model in which
the CO2 and brine fluid properties are now a function of the aquifer conditions,
such as pressure, temperature, and salinity. Specifically:

- The CO2 component is present in the CO2-rich phase but can also dissolve in the brine phase. The amount of dissolved CO2 depends on pressure, temperature, and salinity. For now, in GEOS, the water component cannot be present in the CO2-rich phase.
- Densities and viscosities depend nonlinearly on pressure, temperature, and salinity.
- The hydrostatic initial condition accounts for the geothermal gradient of 0.03 K/m specified in the benchmark description.

We plan to use two types of physical models in this benchmark:

- A model simulating flow and mass transfer, but not heat transfer (i.e., no energy balance is used). The geothermal gradient is constant in time, and is taken into account in the calculation of temperature-dependent properties.
- A fully thermal model simulating flow as well as mass and heat transfer. The results obtained with this more complex model are not available yet and will be added to this page later.

------------------------------------------------------------------
Mesh and element regions
------------------------------------------------------------------

As illustrated by the ECLIPSE results in `(Class et al., 2009) <https://link.springer.com/article/10.1007/s10596-009-9146-x>`__,
the leakage rate exhibits a high dependence on the degree of spatial refinement (particularly between the two wells in our observations).
Therefore, we consider two meshes in this test case:

- A "coarse" mesh with 206070 cells, whose spatial resolution is similar to that used by most codes based on the information provided by Table 13 of `(Class et al., 2009) <https://link.springer.com/article/10.1007/s10596-009-9146-x>`__. 
- A "fine" mesh with 339390 cells, whose spatial resolution is finer between the two wells.   

These structured meshes are defined as in :ref:`ExampleIsothermalLeakyWell`, as shown next for the "fine" mesh.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/thermalLeakyWell/thermalLeakyWell_benchmark.xml
    :language: xml
    :start-after: <!-- SPHINX_MESH -->
    :end-before: <!-- SPHINX_MESH_END -->

As in the previous benchmark, we define four element regions whose material list now includes the name of
the capillary pressure constitutive model (``cappres``).
We refer the reader to :ref:`ExampleIsothermalLeakyWell` for an example of this procedure.

------------------------------------------------------------------
Flow solver
------------------------------------------------------------------

Although the fluid behavior is significantly different from that of the previous benchmark, we still use the GEOS
general-purpose multiphase flow solver defined in the XML block **CompositionalMultiphaseFVM**:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/thermalLeakyWell/thermalLeakyWell_base_iterative.xml
    :language: xml
    :start-after: <!-- SPHINX_SOLVER -->
    :end-before: <!-- SPHINX_SOLVER_END -->

.. note::
   The attribute ``temperature`` listed above is mandatory, but are overridden by GEOS to impose a non-uniform geothermal gradient along the z-axis, as we will see later. 

------------------------------------------------------------------
Constitutive models
------------------------------------------------------------------

The Brooks-Corey relative permeabilities and capillary pressure are described using tables
constructed from the parameters values provided in the benchmark description, with a wetting-phase
saturation range between 0.2 and 0.95, an entry pressure of 10000 Pa, and a Brooks-Corey parameter of 2.
We refer the reader to the files used in the **TableFunction** listed below for the exact values
that we have used:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/thermalLeakyWell/thermalLeakyWell_base_iterative.xml
    :language: xml
    :start-after: <!-- SPHINX_SCAL -->
    :end-before: <!-- SPHINX_SCAL_END -->

The two-phase, two-component CO2-brine model implemented in GEOS is parameterized in
the **CO2BrinePhillips** XML block:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/thermalLeakyWell/thermalLeakyWell_base_iterative.xml
    :language: xml
    :start-after: <!-- SPHINX_FLUID -->
    :end-before: <!-- SPHINX_FLUID_END -->

The components of this fluid model are described in detail in the :ref:`CO2-EOS` and are briefly summarized below.
They are parameterized using three parameter files that must be written carefully
to obtain the desired behavior, as explained next.

CO2 density and viscosity
~~~~~~~~~~~~~~~~~~~~~~~~~

These properties are obtained using the models proposed by Span and Wagner (1996) and Fenghour and Wakeham (1998)
for density and viscosity, respectively.
The density and viscosity values are internally tabulated by GEOS at the beginning of the simulation by solving
the Helmholtz energy equation for each pair :math:`(p,T)`.

The tables size and spacing are specified in the file `pvtgas.txt`.
Here, for both quantities, the values are tabulated between 6.6e6 Pa and 4e7 Pa, with a pressure spacing of 1e6 Pa,
and between 302 K and 312 K, with a temperature increment of 5 K.
These values have been chosen using the initial condition and an upper bound on the expected pressure increase
during the simulation.

.. code:: 

        DensityFun SpanWagnerCO2Density 6.6e6 4e7 1e6 302.0 312.0 5
        ViscosityFun FenghourCO2Viscosity 6.6e6 4e7 1e6 302.0 312.0 5

.. note::
   If pressure or temperature go outside the values specified in this parameter file, constant extrapolation is used to obtain the density and viscosity values. Note that for now, no warning is issued by GEOS when this happens. We plan to add a warning message to document this behavior in the near future.  
	
Brine density and viscosity
~~~~~~~~~~~~~~~~~~~~~~~~~~~

These properties depend on pressure, temperature, composition, and salinity via the models proposed by
Phillips et al. (1981). The brine density is modified to account for the presence of dissolved CO2 using
the method proposed by Garcia (2001).   
The values of (pure) brine density are also tabulated at a function of pressure and temperature, and we
use the same range as for the CO2 properties to construct this table:

.. code::

        DensityFun PhillipsBrineDensity 6.6e6 4e7 1e6 302.0 312.0 5 1.901285269
        ViscosityFun PhillipsBrineViscosity 1.901285269

Importantly, the last value on each line in the file `pvtliquid.txt` defines the salinity in the domain.
In our model, salinity is constant in space and in time (i.e., unlike water and CO2, it is not tracked as
a component in GEOS).
In our model, salinity is specified as a molal concentration in mole of NaCl per kg of solvent (brine).
The value used here (1000 x 10 / ( 58.44 x ( 100 - 10 ) ) = 1.901285269 moles/kg) is
chosen to match the value specified in the benchmark (weight% of 10%).

CO2 solubility in brine
~~~~~~~~~~~~~~~~~~~~~~~

As explained in :ref:`CO2-EOS`, we use the highly nonlinear model proposed by Duan and Sun (2004)
to compute the CO2 solubility as a function of pressure, temperature, composition, and salinity.
In `co2flash.txt`, we use the same parameters as above to construct the pressure-temperature tables
of precomputed CO2 solubility in brine.

.. code::

        FlashModel CO2Solubility 6.6e6 4e7 1e6 302.0 312.0 5 1.901285269
	
------------------------------------------------------------------
Initial and boundary conditions
------------------------------------------------------------------

The domain is initially saturated with brine with a hydrostatic pressure field and a geothermal
gradient of 0.03 K/m.
This is specified using the **HydrostaticEquilibrium** XML tag in the **FieldSpecifications** block:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/thermalLeakyWell/thermalLeakyWell_base_iterative.xml
    :language: xml
    :start-after: <!-- SPHINX_HYDROSTATIC -->
    :end-before: <!-- SPHINX_HYDROSTATIC_END -->

Although this is the same block as in :ref:`ExampleIsothermalLeakyWell`, GEOS is now enforcing the
geothermal gradient specified in the **TableFunction** named ``initTempTable``, and is also accounting
for the nonlinear temperature dependence of brine density to equilibrate the pressure field.

We use the simple table-based approach shown below to impose the Dirichlet boundary conditions on the four
sides of the domain.  

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/thermalLeakyWell/thermalLeakyWell_base_iterative.xml
    :language: xml
    :start-after: <!-- SPHINX_DIRICHLET_BC -->
    :end-before: <!-- SPHINX_DIRICHLET_BC_END -->

where the ``setNames = "{ east, west, south, north }"`` are defined using the **Box**
XML tags of the **Geometry** section, and where the tables are defined as **TableFunction**
in the **Functions** section. 
		 
.. note::
   Due to the nonlinear dependence of brine density on temperature, this block does not exactly impose a Dirichlet pressure equal to the initial condition. Instead, here, we impose a linear pressure gradient along the z-axis, whose minimum and maximum values are the same as in the initial state. We could have imposed Dirichlet boundary conditions preserving the initial condition using as many points in `zlin.geos` as there are cells along the z-axis (instead of just two points).  

The **SourceFlux** is the same as in the previous benchmark case (see :ref:`ExampleIsothermalLeakyWell`).

---------------------------------
Inspecting results
---------------------------------

We request VTK-format output files and use Paraview to visualize the results.
The following figures show the distribution of CO2 saturation and pressure along the slice defined by x = 0 at t = 1,000 days.  

.. _thermalLeakyWell_CO2saturation:
.. figure:: co2_saturation.png
   :align: center
   :width: 500
   :figclass: align-center

   CO2 saturation after 1,000 days

.. _thermalLeakyWell_pressure:
.. figure:: pressure.png
   :align: center
   :width: 500
   :figclass: align-center

   Pressure after 1,000 days

To validate the GEOS results, we consider the metrics used in
`(Class et al., 2009) <https://link.springer.com/article/10.1007/s10596-009-9146-x>`__ as
previously done in :ref:`ExampleIsothermalLeakyWell`.

First, we consider the arrival time of the CO2 plume at the leaky well.
As in `(Class et al., 2009) <https://link.springer.com/article/10.1007/s10596-009-9146-x>`__,
we use the leakage rate threshold of 0.005% to detect the arrival time.
In our numerical tests, the arrival time is highly dependent on the degree of spatial refinement
in the vicinity of the wells and on the time step size, but these parameters are not documented in 
`(Class et al., 2009) <https://link.springer.com/article/10.1007/s10596-009-9146-x>`__.
The next table reports the GEOS arrival time at the leaky well and compares it with the values published in  
`(Class et al., 2009) <https://link.springer.com/article/10.1007/s10596-009-9146-x>`__.

+-----------------------------+---------------------+
| Code                        | Arrival             |
|                             | time [day]          |
+=============================+=====================+
| GEOSX COARSE                |     36.8            | 
+-----------------------------+---------------------+
| GEOSX FINE                  |     46.7            | 
+-----------------------------+---------------------+
| COORES                      |     31              | 
+-----------------------------+---------------------+
| ECLIPSE HW                  |     42              | 
+-----------------------------+---------------------+
| ECLIPSE SCHLUMBERGER COARSE |     24              |
+-----------------------------+---------------------+
| ECLIPSE SCHLUMBERGER FINE   |     34              |
+-----------------------------+---------------------+
| RockFlow                    |     30              | 
+-----------------------------+---------------------+
| TOUGH2                      |     46              | 
+-----------------------------+---------------------+

.. note::
   In the table above, we only included the values obtained with the codes that do **not** solve an energy balance equation. The values obtained with the fully thermal codes (FEHM, MUFTE, and RTAFF2) are omitted for now.

Next, we measure the CO2 leakage rate through the leaky well, defined by the authors as the CO2 mass
flow at midway between top and bottom aquifers divided by the injection rate (8.87 kg/s), in percent.
The GEOS leakage rate is shown in the figure below:

.. plot:: docs/sphinx/advancedExamples/validationStudies/carbonStorage/thermalLeakyWell/thermalLeakyWell.py

We see that GEOS produces a reasonable match with the numerical codes considered in the study.
Although it is not possible to exactly match the published results (due to the lack of information
on the problem, such as mesh refinement and time step size), GEOS reproduces well the trend exhibited
by the other codes.
	
For reference, we include below the original figure from
`(Class et al., 2009) <https://link.springer.com/article/10.1007/s10596-009-9146-x>`__
containing all the results, including those obtained with the codes solving an energy equation.

.. _thermalLeakyWell_referenceLeakageRate:
.. figure:: referenceLeakageRates.png
   :align: center
   :width: 500
   :figclass: align-center

   Leakage rates [%] obtained with the simulators considered in `(Class et al., 2009) <https://link.springer.com/article/10.1007/s10596-009-9146-x>`__.

To further validate the GEOS results, we reproduce below Table 9 of
`(Class et al., 2009) <https://link.springer.com/article/10.1007/s10596-009-9146-x>`__ (only considering
codes that do not solve an energy equation) to compare the maximum leakage rate, the time at which this
maximum leakage rate is attained, and the leakage rate at 2000 days.
We observe that the GEOS values are in the same range as those considered in the benchmark.

+-----------------------------+-------------------------+---------------------------+--------------------------+
| Code                        | Max                     | Time at                   | Leakage at               |
|                             | leakage [%]             | max leakage [day]         | 2000 days [%]            |
+=============================+=========================+===========================+==========================+
| GEOSX COARSE                |  0.102                  |    438.5                  |    0.075                 |
+-----------------------------+-------------------------+---------------------------+--------------------------+
| GEOSX FINE                  |  0.115                  |    425.0                  |    0.085                 |
+-----------------------------+-------------------------+---------------------------+--------------------------+
| COORES                      |  0.105                  |    300                    |    0.076                 |
+-----------------------------+-------------------------+---------------------------+--------------------------+
| ECLIPSE HW                  |  0.074                  |    600                    |    0.067                 |
+-----------------------------+-------------------------+---------------------------+--------------------------+
| ECLIPSE SCHLUMBERGER COARSE |  0.109                  |    437                    |    0.086                 |
+-----------------------------+-------------------------+---------------------------+--------------------------+
| ECLIPSE SCHLUMBERGER FINE   |  0.123                  |    465                    |    0.094                 |
+-----------------------------+-------------------------+---------------------------+--------------------------+
| RockFlow                    |  0.11                   |    279                    |    0.09                  |
+-----------------------------+-------------------------+---------------------------+--------------------------+
| TOUGH2                      |  0.096                  |    400                    |    0.067                 |
+-----------------------------+-------------------------+---------------------------+--------------------------+

This table confirms the agreement between GEOS and the results of `(Class et al., 2009) <https://link.springer.com/article/10.1007/s10596-009-9146-x>`__.

------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
