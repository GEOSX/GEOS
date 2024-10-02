.. _EquilibriumInitialCondition:

####################################################
Hydrostatic Equilibrium Initial Condition
####################################################

Overview
======================

The user can request an initialization procedure enforcing a hydrostatic equilibrium for flow simulations and for coupled flow and mechanics simulations.
The hydrostatic initialization is done by placing one or more **HydrostaticEquilibrium** tag(s) in the **FieldSpecifications** block of the XML input file.
This initialization procedure is described below in the context of single-phase and compositional multiphase flow.
At the end of this document, we compare the hydrostatic equilibrium method to another initialization method, based on the input of x-y-z tables.

Single-phase flow parameters
==============================

For single-phase flow, the **HydrostaticEquilibrium** initialization procedure requires the following user input parameters:

* ``datumElevation``: the elevation (in meters) at which the datum pressure is enforced. The user must ensure that the datum elevation is within the elevation range defined by the input mesh. GEOS issues a warning if this is not the case.

* ``datumPressure``: the pressure value (in Pascal) enforced by GEOS at the datum elevation. 

* ``objectPath``: the path defining the groups on which the hydrostatic equilibrium is computed. We recommend using ``ElementRegions`` to apply the hydrostatic equilibrium to all the cells in the mesh. Alternatively, the format ``ElementRegions/NameOfRegion/NameOfCellBlock`` can be used to select only a cell block on which the hydrostatic equilibrium is computed.

.. note::
   In GEOS, the z-axis is positive going upward, this is why the attributes listed in this page are expressed as a function of elevation, not depth. 

Using these parameters and the pressure-density constitutive relationship, GEOS uses a fixed-point iteration scheme to populate a table of hydrostatic pressures as a function of elevation. The fixed-point iteration scheme uses two optional attributes: ``equilibriumTolerance``, the absolute tolerance to declare that the algorithm has converged, and ``maxNumberOfEquilibrationTolerance``, the maximum number of iterations for a given elevation in the fixed point iteration scheme.

In addition, the elevation spacing of the hydrostatic pressure table is set with the optional ``elevationIncrementInHydrostaticPressureTable`` parameter (in meters), whose default value is 0.6096 meters. 
Then, once the table is fully constructed, the hydrostatic pressure in each cell is obtained by interpolating in the hydrostatic pressure table using the elevation at the center of the cell.

.. note::
   The initialization algorithm assumes that the ``gravityVector`` (defined in the **Solvers** XML tag) is aligned with the z-axis. If this is not the case, GEOS terminates the simulation when the **HydrostaticEquilibrium** tag is detected in the XML file. 

Compositional multiphase flow parameters
==========================================

For compositional multiphase flow, the **HydrostaticEquilibrium** initialization procedure follows the same logic but requires more input parameters.
In addition to the required ``datumElevation``, ``datumPressure``, and ``objectPath`` parameters listed above, the user must specify:

* ``componentNames``: the names of the components present in the fluid model. This field is used to make sure that the components provided to **HydrostaticEquilibrium** are consistent with the components listed in the fluid model of the **Constitutive** block. 

* ``componentFractionVsElevationTableNames``: the names of :math:`n_c` tables (where :math:`n_c` is the number of components) specifying the component fractions as a function of elevation. There must be one table name per component, and the table names must be listed in the same order as the components in ``componentNames``. 

* ``temperatureVsElevationTableName``: the names of the table specifying the temperature (in Kelvin) as a function of elevation.

* ``initialPhaseName``: the name of the phase initially saturating the domain. The other phases are assumed to be at residual saturation at the beginning of the simulation. 

These parameters are used with the fluid density model (depending for compositional flow on pressure, component fractions, and in some cases, temperature) to populate the hydrostatic pressure table, and later initialize the pressure in each cell.

.. note::
   The current initialization algorithm has an important limitation and does not support initial phase contacts (e.g., water-oil, gas-oil, or water-gas contacts). The implementation assumes only one mobile phase in the initial system, identified by the ``initialPhaseName`` attribute. The other phases are assumed at residual saturation. As a result, the system may not be at equilibrium if there is initially more than one mobile phase in the system (for instance if the domain is saturated with gas at the top, and water at the bottom, for instance).

.. note::
   As in the single-phase flow case, GEOS terminates the simulation if **HydrostaticEquilibrium** tag is present in an XML file defining a ``gravityVector`` not aligned with the z-axis.

The full list of parameters is provided below:

.. include:: /docs/sphinx/datastructure/HydrostaticEquilibrium.rst


Examples
=======================

For single-phase flow, a typical hydrostatic equilibrium input looks like:

.. code-block:: xml
	     
   <FieldSpecifications>
   
      <HydrostaticEquilibrium
        name="equil"
        objectPath="ElementRegions"      
        datumElevation="5"
        datumPressure="1e6"/>
      
   </FieldSpecifications>

For compositional multiphase flow, using for instance the CO2-brine flow model, a typical hydrostatic equilibrium input looks like:

.. code-block:: xml

   <FieldSpecifications>		
	     
      <HydrostaticEquilibrium
        name="equil"
        objectPath="ElementRegions"      
        datumElevation="28.5"
        datumPressure="1.1e7"
        initialPhaseName="water"
        componentNames="{ co2, water }"
        componentFractionVsElevationTableNames="{ initCO2CompFracTable,
                                                  initWaterCompFracTable }"
        temperatureVsElevationTableName="initTempTable"/>

   </FieldSpecifications>

In this case, a possible way to provide the three required tables is:

.. code-block:: xml

   <Functions>

     <TableFunction
       name="initCO2CompFracTable"
       coordinates="{ 0.0, 10.0, 20.0, 30.0 }"
       values="{ 0.04, 0.045, 0.05, 0.055 }"/>

     <TableFunction
       name="initWaterCompFracTable"
       coordinates="{ 0.0, 10.0, 20.0, 30.0 }"
       values="{ 0.96, 0.955, 0.95, 0.945 }"/>

     <TableFunction
       name="initTempTable"
       coordinates="{ 0.0, 15.0, 30.0 }"
       values="{ 358.15, 339.3, 333.03 }"/>
     
   </Functions>

Note that the spacing of the two component fraction tables must be the same, but the spacing of the temperature table can be different.

Expected behavior and comparison with another initialization method
=====================================================================

As illustrated in :ref:`TutorialFieldCase`, users can also use multiple **FieldSpecification** tags to impose initial fields, such as the pressure, component fractions, and temperature fields.
To help users select the initialization method that best meets their needs, we summarize and compare below the two possible ways to initialize complex, non-uniform initial fields for compositional multiphase simulations in GEOS.

Initialization using **HydrostaticEquilibrium**
-----------------------------------------------

This is the initialization procedure that we have described in the first sections of this page.
In **HydrostaticEquilibrium**, the initial component fractions and temperatures are provided as a function of elevation only, and the hydrostatic pressure is computed internally before the simulation starts.
The typical input was illustrated for a CO2-brine fluid model in the previous paragraph.

Expected behavior:

* If **FieldSpecification** tags specifying initial pressure, component fractions, and/or temperature are included in an XML input file that also contains the **HydrostaticEquilibrium** tag, the **FieldSpecification** tags are ignored by GEOS. In other words, only the pressure, component fractions, and temperature fields defined with the **HydrostaticEquilibrium** tag as a function of elevation are taken into account.

* In the absence of source/sink terms and wells, the initial flow residual should be smaller than :math:`10^-6`. Similarly, in coupled simulations, the residual of the mechanical problem should be close to zero.

Initialization using **FieldSpecification** tags
------------------------------------------------

This is the initialization method illustrated in :ref:`TutorialFieldCase`.
The user can impose initial pressure, component fractions, and temperature fields using **FieldSpecification** tags, such as, for a two-component CO2-brine case:

.. code-block:: xml

   <FieldSpecifications>		
	     
     <FieldSpecification
       name="initialPressure"
       initialCondition="1"
       setNames="{ all }"
       objectPath="ElementRegions"
       fieldName="pressure"
       scale="1"
       functionName="initialPressureTableXYZ"/>

     <FieldSpecification
       name="initialCO2CompFraction"
       initialCondition="1"
       setNames="{ all }"
       objectPath="ElementRegions"
       fieldName="globalCompFraction"
       component="0"
       scale="1"
       functionName="initialCO2CompFracTableXYZ"/>

     <FieldSpecification
       name="initialWaterCompFrac"
       initialCondition="1"
       setNames="{ all }"
       objectPath="ElementRegions"
       fieldName="globalCompFraction"
       component="1"
       scale="1"
       functionName="initialWaterCompFracTableXYZ"/>

     <FieldSpecification
       name="initialTemperature"
       initialCondition="1"
       setNames="{ all }"
       objectPath="ElementRegions"
       fieldName="temperature"
       scale="1"
       functionName="initialTemperatureTableXYZ"/>
       
   </FieldSpecifications>		       

In this input method, ``initialPressureTableXYZ``, ``initialCO2CompFracTableXYZ``, ``initialWaterCompFracTableXYZ``, and ``initialTemperatureTableXYZ`` are tables describing these initial fields as a function of the x, y, and z spatial coordinates.
Then, the cell-wise values are determined by interpolating in these tables using the coordinates of the center of the cells.

Expected behavior:

* In this approach, it is the responsibility of the user to make sure that these initial fields satisfy a hydrostatic equilibrium. If not, the model will equilibrate itself during the first time steps of the simulation, possibly causing large initial changes in pressure, component fractions, and temperature.

* If the initial state imposed by the **FieldSpecification** tags is not at equilibrium, the displacements produced by coupled flow and mechanics simulations should be interpreted with caution, as these displacements are computed with respect to a non-equilibrium initial state.

* This method is suited to impose initial fields in complex cases currently not supported by the **HydrostaticEquilibrium** tag (e.g., in the presence of phase contacts, capillary pressure, etc). Specifically, the user can equilibrate the model using other means (such as using another simulator, or running a few steps of GEOS), retrieve the equilibrated values, convert them into x-y-z tables, and impose them in the new GEOS simulations using **FieldSpecification** tags.  
