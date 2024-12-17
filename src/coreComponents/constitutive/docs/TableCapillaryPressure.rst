.. _TableCapillaryPressure:

############################################
Table capillary pressure
############################################


Overview
======================

The user can specify the capillary pressures using tables describing a piecewise-linear capillary pressure function of volume fraction (i.e., saturation) for each phase, except the reference phase for which capillary pressure is assumed to be zero.
Depending on the number of fluid phases, this model is used as follows:

* For two-phase flow, the user must specify one capillary pressure table. During the simulation, the capillary pressure of the non-reference phase is computed by interpolating in the table as a function of the non-reference phase saturation.   

* For three-phase flow, the user must specify two capillary pressure tables. One capillary pressure table is required for the pair wetting-phase--intermediate-phase (typically, water-oil), and one capillary pressure table is required for the pair non-wetting-phase--intermediate-phase (typically, gas-oil). During the simulation, the former is used to compute the wetting-phase capillary pressure as a function of the wetting-phase volume fraction and the latter is used to compute the non-wetting-phase capillary pressure as a function of the non-wetting-phase volume fraction. The intermediate phase is assumed to be the reference phase, and its capillary pressure is set to zero.    

Below is a table summarizing the choice of reference pressure for the various phase combinations:
  
============================ ================
Phases present in the model  Reference phase
============================ ================
oil, water, gas              Oil phase 
oil, water                   Oil phase
oil, gas                     Oil phase
water, gas                   Gas phase
============================ ================

In all cases, the user-provided capillary pressure is used in GEOS to compute the phase pressure using the formula:

.. math::
    P_c = p_{nw} - p_w.

where :math:`p_{nw}` and :math:`p_w` are respectively the non-wetting-phase and wetting-phase pressures.    
  
Parameters
======================

The capillary pressure constitutive model is listed in
the ``<Constitutive>`` block of the input XML file.
The capillary pressure model must be assigned a unique name via
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute of the ``<ElementRegions>``
node.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/TableCapillaryPressure.rst

Below are some comments on the model parameters.

* ``phaseNames`` - The number of phases can be either two or three. Supported phase names are:

===== ===========
Value Phase
===== ===========
oil   Oil phase
gas   Gas phase
water Water phase
===== ===========

* ``wettingNonWettingCapPressureTableName`` - The name of the capillary pressure table for two-phase systems. Note that this keyword is only valid for two-phase systems, and is not allowed for three-phase systems (for which the user must specify instead ``wettingIntermediateCapPressureTableName`` and ``nonWettingIntermediateCapPressureTableName``). This capillary pressure must be a strictly decreasing function of the water-phase volume fraction (for oil-water systems and gas-water systems), or a strictly increasing function of the gas-phase volume fraction (for oil-gas systems).    

* ``wettingIntermediateCapPressureTableName`` - The name of the capillary pressure table for the pair wetting-phase--intermediate-phase. This capillary pressure is applied to the wetting phase, as a function of the wetting-phase volume fraction. Note that this keyword is only valid for three-phase systems, and is not allowed for two-phase systems (for which the user must specify instead ``wettingNonWettingCapPressureTableName``). This capillary pressure must be a strictly decreasing function of the wetting-phase volume fraction.

* ``nonWettingIntermediateCapPressureTableName`` - The name of the capillary pressure table for the pair non-wetting-phase--intermediate-phase. Note that this keyword is only valid for three-phase systems, and is not allowed for two-phase systems (for which the user must specify instead ``wettingNonWettingCapPressureTableName``). This capillary pressure must be a strictly increasing function of the non-wetting-phase volume fraction. 


Examples
=======================

For a two-phase water-gas system (for instance in the CO2-brine fluid model), a typical capillary pressure input looks like:

.. code-block:: xml

   <Constitutive>
      ...
      <TableCapillaryPressure
        name="capPressure"
        phaseNames="{ water, gas }"
        wettingNonWettingCapPressureTableNames="waterCapillaryPressureTable"/>
      ...
   </Constitutive>

For a three-phase oil-water-gas system (for instance in the Black-Oil fluid model), a typical capillary pressure input looks like:

.. code-block:: xml

   <Constitutive>
      ...
      <TableCapillaryPressure
        name="capPressure"
        phaseNames="{ water, oil, gas }"
        wettingIntermediateCapPressureTableName="waterCapillaryPressureTable"
        nonWettingIntermediateCapPressureTableName="gasCapillaryPressureTable"/>	
      ...
   </Constitutive>

The tables mentioned above by name must be defined in the ``<Functions>`` block of the XML file using the ``<TableFunction>`` keyword. 
