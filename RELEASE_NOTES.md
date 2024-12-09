[comment]: # (-----------------------------------------------------------------)
[comment]: # (SPDX-License-Identifier: LGPL-2.1-only)
[comment]: # 
[comment]: # (Copyright 2018-2020 Lawrence Livermore National Security LLC)
[comment]: # (Copyright 2018-2020 The Board of Trustees of the Leland Stanford)
[comment]: # (                    Junior University)
[comment]: # (Copyright 2018-2020 TotalEnergies)
[comment]: # (Copyright 2019-     GEOS/GEOSX Contributors)
[comment]: # (All right reserved)
[comment]: # 
[comment]: # (For more details see:)
[comment]: # (  https://github.com/GEOS-DEV/GEOS/LICENSE)
[comment]: # (  https://github.com/GEOS-DEV/GEOS/COPYRIGHT)
[comment]: # (  https://github.com/GEOS-DEV/GEOS/CONTRIBUTORS)
[comment]: # (  https://github.com/GEOS-DEV/GEOS/NOTICE)
[comment]: # (  https://github.com/GEOS-DEV/GEOS/ACKNOWLEDGEMENTS)
[comment]: # (  https://github.com/GEOS-DEV/GEOS/RELEASE)


Version v0.2.0 -- Release date 2020-06-20
==========================================
* Extensive updates to Data Repository
* Discretization
  * Finite element interface
    * Kernel launching interface looping abstraction
    * Element formulations for 1st order:
      * 8-node hexahedron
      * 6-node wedge
      * 5-node pyramid
      * 4-node tetrahedron
  * Cell-centered Finite Volume method with Two-point Flux Approximation (TPFA)
  * Hybrid Finite Volume method with TPFA and quasi-TPFA inner products
* Physics Solvers
  * Solid mechanics explicit on GPU, implicit assembly on GPU
  * Single-phase flow (assembly on GPU)
    * Classical FVM and Hybrid FVM formulations
    * Porous matrix and DFM fracture flow
  * Compositional multiphase flow (assembly on GPU)
    * Fully implicit isothermal overall composition formulation
    * Fluid constitutive models:
      * Equation-of-state hydrocarbon compositional
      * Three-phase extended black-oil
      * Two-phase CO2-brine
  * Multi-segmented wells for single phase and compositional multiphase flow (assembly on GPU)
  * Surface Generation
    * Topology change (legacy GEOS approach)
* Embedded Discrete Fractures
  * Enriched finite element method for the discretization of the mechanics
  * Piecewise constant displacement jump enrichment
  * Hydrofracture solver (legacy GEOS approach)
  * Small strain aligned contact using Lagrange multipliers
    * Discrete fracture model using a low-order stabilized mixed finite element method
  * Proppant Transport Solver
    * FVM formulation
    * Major physical processes modeled:
      * Proppant-fluid slurry flow and multicomponent transport in fractures
      * Proppant gravitational settling
      * Proppant bed build-up and development space
* Mesh Structure
  * Introduced the concept of extrinsic mesh data
  * Fracture elements to represent FV cells in fractures
* VTK output
* Linear algebra interface layers for Hypre, Trilinos, Petsc
  * Common interface for supported linear algebra packages
  * Krylov solvers (CG, GMRES, BiCGSTAB)
  * Preconditioners (algebraic multigrid, incomplete factorizations)
  * Block matrix and vector support
  * Serial and parallel direct solvers


Version v0.1.0 -- Release date 2018-02-15
==========================================
Initial Code Release containing:
* Data Repository
  * Group
  * Wrapper
  * Input processing 
* Physics solver hierarchy
  * Solid Mechanics
* Mesh data structure
  * NodeManager, EdgeManager, FaceManager, ElementManager
  * Silo Output
