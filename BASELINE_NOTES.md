
Notes
==========

This file is designed to track changes to the integrated test baselines.
Any developer who updates the baseline ID in the .integrated_tests.yaml file is expected to create an entry in this file with the pull request number, date, and their justification for rebaselining.
These notes should be in reverse-chronological order, and use the following time format: (YYYY-MM-DD).

PR #3339 (2024-11-08)
====================
Fix default parameter (amgNumFunctions = 3) for models involving mechanics.

PR #3337 (2024-11-06)
====================
Change density treatment for gravity in multiphase flow solver.

PR #3408 (2024-11-06)
====================
EFEM bugfixes: effective traction + oldStress.

PR #3280 (2024-11-05)
====================
Added Sprig-slider test.

PR #2909 (2024-10-30)
=====================
Add routine for automatic time steps in waveSolvers with new attributes

PR #3156 (2024-10-29)
====================
Restart check errors due to 1) schema node added to enable thermal option in well model and 2) arrays removed/added for option.  Max difference errors due treatment of shutin wells.  Previously non-zero rate value reported for shutin well, new code will set rate arrays to zero. 

PR #2878 (2024-10-17)
=====================
Sorted region cellBlocks names alphabetically. Therefore affected ordering of: faceManager/elemSubRegionList, nodeManager/elemList, nodeManager/elemSubRegionList, SurfaceElementSubRegion::fractureElementsToCellSubRegions, field::perforation::reservoirElementSubregion.

PR #3364( 2024-10-15)
=====================
Enable reservoir+wells+contact mechanics. Rebaseline needed because of 'allowNegativePressure' flag added for wells.

PR #3364( 2024-10-01)
=====================
Separate mass and volume residuals for output in compositional flow solver. Baseline update because of minor numerical diffs.

PR #3149( 2024-09-30)
=====================
Added new field "writeCSV"

PR #3163 (2024-09-20)
=====================
Added new fields (krylovStrongestTol, adaptiveGamma, adaptiveExponent) to the LinearSolverParameters for adaptive tolerances.

PR #3338 (2024-09-19)
======================
Updated time-stepping logic. Rebaseline due to new input parameter and minor numerical diffs.

PR #3217 (2024-09-16)
======================
ALM slip and open modes with relative tests.

PR #3318 (2024-09-12)
======================
Modified SeismicityRate poroelastic case.

PR #3322 (2024-09-06)
======================
Print out fracture state for contact model. Rebaseline the corresponding cases.

PR #3302 (2024-09-05)
======================
Added restartcheks to hydrofrac cases and reduced time of cases that were too long.

PR #3135 (2024-09-04)
======================
Temperature dependent single phase thermal conductivity. Rebaseline all thermal cases.

PR #3294 (2024-09-01)
======================
Re-enable enforcement of wave propagation integrated test pass.

PR #3300 (2024-08-28)
======================
Re-enable floating point exceptions. Rebaseline due to minor changing default value of maxRelativeCompDensChange from 1.7976931348623157e+308 to 1.7976931348623157e+208.

PR #3283 (2024-08-22)
======================
Reuse computeSinglePhaseFlux. Rebaseline due to minor numerical diffs.

PR #3249 (2024-08-14)
======================
Two initialization options for poromechanical models. Rebaseline the corresponding cases.

PR #3278 (2024-08-12)
======================
Renamed GEOSX to GEOS in enternal mesh import, so rebaseline to change these names is the baselines.

202 (2024-08-03)
======================
Acoustic VTI tests needed rebaselining after update in source and receiver location algorithm.

PR #3215 (2024-07-23)
======================
Changed the default value for massCreation and name of the wrapper.

PR #3194 (2024-07-22)
======================
Check pore volume for all element types, also check that default aperture > 0. Rebaseline for modified tests. No real results change.

PR #3213 (2024-07-12)
======================
Added baselines for new tests on Dirichlet boundary conditions for multiphase flow.

PR #3194 (2024-07-10)
======================
Use aperture table in poromechanics with conforming fractures. Rebaseline the corresponding cases.

PR #3006 (2024-07-01)
======================
Added baselines for new tests. Relaxing tolerances for singlePhasePoromechanics_FaultModel_smoke.

PR #3196 (2024-06-28)
======================
Added isLaggingFractureStencilWeightsUpdate to hydrofracture solve. Rebaseline because of the new input.

PR #3177 (2024-06-28)
======================
Added logLevel to TimeHistoryOutput. Rebaseline because of the new input flag.

PR #3181 (2024-06-25)
======================
Decouple debug matrix output from logLevel. Rebaseline because of the new input flag.

PR #3142 (2024-06-20)
======================
Adding output of total strain. Rebaseline because of new inclusion of strain in output.

PR #3170 (2024-06-19)
======================
Fix tutorial example for thermal debonding wellbore problem. Test case modified.

PR #3130 (2024-06-19)
======================
New solver for contact mechanics based on the Augmented Lagrangian Method (ALM). New test case added.

PR #3160 (2024-06-18)
======================
Two experimental options for compositional flow solver. Rebaseline because of the new input flags.

PR #3165 (2024-06-18)
======================
Small bug fix. Rebaseline required due to appearance of useTotalMassEquation in well solver params. No real results change.

PR #3088 (2024-06-17)
======================
Adding temperature-dependent Solid Volumetric Heat Capacity. Rebaseline because of the parameter change in SolidInternalEnergy.

PR #3100 (2024-06-14)
======================
Adding pressure stabilization for single phase poromechanics.

PR #3133 (2024-06-14)
======================
Fix node ordering for faceElements.

PR #3021 (2024-06-13)
======================
Preparatory work for fractures + wells. New test case added.

PR #3152 (2024-06-13)
======================
Some random things. Baseline update because of the new parameter (minScalingFactor).

PR #3138 (2024-06-11)
======================
Properly sync nonlinear solver params for coupled solver. Baseline update mostly due to number of iterations change in baseline files.

PR #3140 (2024-06-11)
======================
Fixed derivative in EzrokhiBrineDensity

PR #3080 (2024-06-07)
=====================
Rebaseline after adding viscoelastic wave propagator.

PR #3075 (2024-06-05)
=====================
Introduce configuration tolerance. Rebaseline because of the new parameter in NonlinearSolverParameters.

PR #3120 (2024-06-05)
=====================
Add missing compositionalMultiphaseFlow tests into ATS and adjust output naming. Rebaseline accordingly.

PR #3113 (2024-06-05)
=====================
Add general version updateConfiguration. Rebaseline of some edfm cases is needed.

PR #3050 (2024-05-20)
=====================
Spatially varying grain bulk modulus. Rebaseline of all poromechanics cases needed.

PR #3141 (2024-05-28)
=====================
Test cashing baselines locally.

PR #3125 (2024-05-16)
=====================
Remove field to store pressure gradient cell-wise for solvers that don't need it.

PR #2110 (2024-05-13)
=====================
new field to store pressure gradient cell-wise.

PR #3060 (2024-05-13)
======================
Rebaselined after addition of elastic VTI wave propagator.

PR #4950 (2024-05-10)
======================
Added smoke tests for SeismicityRate solver in inducedSeismicity.

PR #3086 (2024-05-09)
======================
Added a presure-dependent permeability model and the transmissibility calculation in the CellElementStencil

PR #3105 (2024-05-08)
======================
Added missing derivative for temperature, hence small numerical diffs in thermal tests results and numeracal behavior

PR #2917 (2024-05-07)
======================
New fields for wellsControls: wellControls1_ConstantMassRate_table, targetMassRate, massDensity, ...

PR #3044 (2024-05-02)
======================
Removed old integratedTests submodule
Implemented new baseline storage
Implemented new CI integrated tests
