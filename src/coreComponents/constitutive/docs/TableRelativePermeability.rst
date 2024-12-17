.. _TableRelativePermeability:

############################################
Table relative permeability 
############################################


Overview
======================

The user can specify the relative permeabilities using tables describing a piecewise-linear relative permeability function of volume fraction (i.e., saturation) for each phase.
Depending on the number of fluid phases, this model is used as follows:

* For two-phase flow, the user must specify two relative permeability tables, that is, one for the wetting-phase relative permeability, and one for the non-wetting phase relative permeability. During the simulation, the relative permeabilities are then obtained by interpolating in the tables as a function of phase volume fraction. 

* For three-phase flow, following standard reservoir simulation practice, the user must specify four relative permeability tables. Specifically, two relative permeability tables are required for the pair wetting-phase--intermediate phase (typically, water-oil), and two relative permeability tables are required for the pair non-wetting-phase--intermediate phase (typically, gas-oil). During the simulation, the relative permeabilities of the wetting and non-wetting phases are computed by interpolating in the tables as a function of their own phase volume fraction. The intermediate phase relative permeability is obtained by interpolating the two-phase relative permeabilities using the Baker interpolation procedure.  
  
Parameters
======================

The relative permeability constitutive model is listed in
the ``<Constitutive>`` block of the input XML file.
The relative permeability model must be assigned a unique name via
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute of the ``<ElementRegions>``
node.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/TableRelativePermeability.rst

Below are some comments on the model parameters.

* ``phaseNames`` - The number of phases can be either two or three. For three-phase flow, this model applies a Baker interpolation to the intermediate phase relative permeability. Supported phase names are:

===== ===========
Value Phase
===== ===========
oil   Oil phase
gas   Gas phase
water Water phase
===== ===========

* ``wettingNonWettingRelPermTableNames`` - The list of relative permeability table names for two-phase systems, starting with the name of the wetting-phase relative permeability table, followed by the name of the non-wetting phase relative permeability table. Note that this keyword is only valid for two-phase systems, and is not allowed for three-phase systems (for which the user must specify instead ``wettingIntermediateRelPermTableNames`` and ``nonWettingIntermediateRelPermTableNames``).  

* ``wettingIntermediateRelPermTableNames`` - The list of relative permeability table names for the pair wetting-phase--intermediate-phase, starting with the name of the wetting-phase relative permeability table, and continuing with the name of the intermediate phase relative permeability table. Note that this keyword is only valid for three-phase systems, and is not allowed for two-phase systems (for which the user must specify instead ``wettingNonWettingRelPermTableNames``).  

* ``nonWettingIntermediateRelPermTableNames`` - The list of relative permeability table names for the pair non-wetting-phase--intermediate-phase, starting with the name of the non-wetting-phase relative permeability table, and continuing with the name of the intermediate phase relative permeability table. Note that this keyword is only valid for three-phase systems, and is not allowed for two-phase systems (for which the user must specify instead ``wettingNonWettingRelPermTableNames``).  

.. note::     
   We remind the user that the relative permeability must be a strictly increasing function of phase volume fraction. GEOS throws an error when this condition is not satisfied.

Examples
=======================

For a two-phase water-gas system (for instance in the CO2-brine fluid model), a typical relative permeability input looks like:

.. code-block:: xml

   <Constitutive>
      ...
      <TableRelativePermeability
        name="relPerm"
        phaseNames="{ water, gas }"
        wettingNonWettingRelPermTableNames="{ waterRelativePermeabilityTable, gasRelativePermeabilityTable }"/>
      ...
   </Constitutive>

.. note::   
   The name of the wetting-phase relative permeability table must be specified before the name of the non-wetting phase relative permeability table.   

For a three-phase oil-water-gas system (for instance in the Black-Oil fluid model), a typical relative permeability input looks like:

.. code-block:: xml

   <Constitutive>
      ...
      <TableRelativePermeability
        name="relPerm"
        phaseNames="{ water, oil, gas }"
        wettingIntermediateRelPermTableNames="{ waterRelativePermeabilityTable, oilRelativePermeabilityTableForWO }"
        nonWettingIntermediateRelPermTableNames="{ gasRelativePermeabilityTable, oilRelativePermeabilityTableForGO }"/>	
      ...
   </Constitutive>

.. note::   
   For the wetting-phase--intermediate-phase pair, the name of the wetting-phase relative permeability table must be specified first. For the non-wetting-phase--intermediate-phase pair, the name of the non-wetting-phase relative permeability table must be specified first. If the results look incoherent, this is something to double-check.   

   
The tables mentioned above by name must be defined in the ``<Functions>`` block of the XML file using the ``<TableFunction>`` keyword. 
