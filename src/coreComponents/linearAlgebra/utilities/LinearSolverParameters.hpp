/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LinearSolverParameters.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_

#include "codingUtilities/EnumStrings.hpp"

namespace geos
{

/**
 * @brief Set of parameters for a linear solver or preconditioner.
 *
 * This class holds a simple tree of linear solver options.
 * They are set to default values, but can be overwritten as needed.
 */
struct LinearSolverParameters
{
  /**
   * @brief Linear solver type.
   */
  enum class SolverType : integer
  {
    direct,        ///< Direct solver
    cg,            ///< CG
    gmres,         ///< GMRES
    fgmres,        ///< Flexible GMRES
    bicgstab,      ///< BiCGStab
    preconditioner ///< Preconditioner only
  };

  /**
   * @brief Preconditioner type.
   */
  enum class PreconditionerType : integer
  {
    none,      ///< No preconditioner
    jacobi,    ///< Jacobi smoothing
    l1jacobi,  ///< l1-Jacobi smoothing
    fgs,       ///< Gauss-Seidel smoothing (forward sweep)
    sgs,       ///< Symmetric Gauss-Seidel smoothing
    l1sgs,     ///< l1-Symmetric Gauss-Seidel smoothing
    chebyshev, ///< Chebyshev polynomial smoothing
    iluk,      ///< Incomplete LU with k-level of fill
    ilut,      ///< Incomplete LU with thresholding
    ic,        ///< Incomplete Cholesky
    ict,       ///< Incomplete Cholesky with thresholding
    amg,       ///< Algebraic Multigrid
    mgr,       ///< Multigrid reduction (Hypre only)
    block,     ///< Block preconditioner
    direct,    ///< Direct solver as preconditioner
    bgs,       ///< Gauss-Seidel smoothing (backward sweep)
  };

  integer logLevel = 0;     ///< Output level [0=none, 1=basic, 2=everything]
  integer dofsPerNode = 1;  ///< Dofs per node (or support location) for non-scalar problems
  bool isSymmetric = false; ///< Whether input matrix is symmetric (may affect choice of scheme)
  integer stopIfError = 1;  ///< Whether to stop the simulation if the linear solver reports an error

  SolverType solverType = SolverType::direct;          ///< Solver type
  PreconditionerType preconditionerType = PreconditionerType::iluk;  ///< Preconditioner type

  /// Direct solver parameters: used for SuperLU_Dist interface through hypre and PETSc
  struct Direct
  {
    /**
     * @brief How to permute the columns
     */
    enum class ColPerm : integer
    {
      none,        ///< natural
      MMD_AtplusA, ///< multiple minimum degree on At+A
      MMD_AtA,     ///< multiple minimum degree on At*A (heavy)
      colAMD,      ///< approximate minimum degree on columns
      metis,       ///< using METIS
      parmetis     ///< using ParMETIS
    };

    /**
     * @brief How to permute the rows
     */
    enum class RowPerm : integer
    {
      none, ///< natural
      mc64  ///< using HSL routine MC64
    };

    integer checkResidual = 0;        ///< Whether to check the linear system solution residual
    integer equilibrate = 1;          ///< Whether to scale the rows and columns of the matrix
    ColPerm colPerm = ColPerm::metis; ///< Columns permutation
    RowPerm rowPerm = RowPerm::mc64;  ///< Rows permutation
    integer replaceTinyPivot = 1;     ///< Whether to replace tiny pivots by sqrt(epsilon)*norm(A)
    integer iterativeRefine = 1;      ///< Whether to perform iterative refinement
    integer parallel = 1;             ///< Whether to use a parallel solver (instead of a serial one)
  }
  direct;                             ///< direct solver parameter struct

  /// Krylov-method parameters
  struct Krylov
  {
    real64 relTolerance = 1e-6;       ///< Relative convergence tolerance for iterative solvers
    integer maxIterations = 200;      ///< Max iterations before declaring convergence failure
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
    integer maxRestart = 50;          ///< Max number of vectors in Krylov basis before restarting (GPUs)
#else
    integer maxRestart = 200;         ///< Max number of vectors in Krylov basis before restarting (CPUs)
#endif
    integer useAdaptiveTol = false;   ///< Use Eisenstat-Walker adaptive tolerance
    real64 weakestTol = 1e-3;         ///< Weakest allowed tolerance when using adaptive method
    real64 strongestTol = 1e-8;       ///< Strongest allowed tolerance when using adaptive method
    real64 adaptiveGamma = 0.1;       ///< Gamma parameter for adaptive method
    real64 adaptiveExponent = 1.0;    ///< Exponent parameter for adaptive method
  }
  krylov;                             ///< Krylov-method parameter struct

  /// Matrix-scaling parameters
  struct Scaling
  {
    integer useRowScaling = false;    ///< Apply row scaling
    integer useRowColScaling = false; ///< Apply row and column scaling (not yet implemented)
  }
  scaling;                            ///< Matrix-scaling parameter struct

  /// Algebraic multigrid parameters
  struct AMG
  {
    /// AMG cycle type
    enum class CycleType : integer
    {
      V,                  ///< V-cycle
      W                   ///< W-cycle
    };

    /// AMG pre/post smoothing option
    enum class PreOrPost : integer
    {
      pre,                ///< pre-smoothing only
      post,               ///< post-smoothing only
      both                ///< pre- and post-smoothing
    };

    /// AMG smoother type
    enum class SmootherType : integer
    {
      default_,           ///< Use LAI's default option
      jacobi,             ///< Jacobi smoothing
      l1jacobi,           ///< l1-Jacobi smoothing
      fgs,                ///< Gauss-Seidel smoothing (forward sweep)
      bgs,                ///< Gauss-Seidel smoothing (backward sweep)
      sgs,                ///< Symmetric Gauss-Seidel smoothing
      l1sgs,              ///< l1-Symmetric Gauss-Seidel smoothing
      chebyshev,          ///< Chebyshev polynomial smoothing
      ilu0,               ///< ILU(0)
      ilut,               ///< Incomplete LU with thresholding
      ic0,                ///< Incomplete Cholesky
      ict                 ///< Incomplete Cholesky with thresholding
    };

    /// AMG coarse solver type
    enum class CoarseType : integer
    {
      default_,           ///< Use LAI's default option
      jacobi,             ///< Jacobi (GPU support in hypre)
      l1jacobi,           ///< l1-Jacobi (GPU support in hypre)
      fgs,                ///< Gauss-Seidel (forward sweep)
      sgs,                ///< Symmetric Gauss-Seidel
      l1sgs,              ///< l1-Symmetric Gauss-Seidel
      chebyshev,          ///< Chebyshev polynomial (GPU support in hypre)
      direct,             ///< Direct solver as preconditioner
      bgs                 ///< Gauss-Seidel smoothing (backward sweep)
    };

    /// AMG coarsening types (HYPRE only)
    enum class CoarseningType : integer
    {
      default_,           ///< Use LAI's default option
      CLJP,               ///< A parallel coarsening algorithm using independent sets
      RugeStueben,        ///< Classical Ruge-Stueben on each processor, followed by a third pass
      Falgout,            ///< Ruge-Stueben followed by CLJP
      PMIS,               ///< Parallel coarsening as CLJP but with lower complexities (GPU support)
      HMIS                ///< Hybrid PMIS coarsening
    };

    /// AMG interpolation type (HYPRE only)
    enum class InterpType : integer
    {
      default_,           ///< Use LAI's default option
      modifiedClassical,  ///< Modified classical
      direct,             ///< Direct (GPU support)
      multipass,          ///< Multipass (GPU support)
      extendedI,          ///< Extended+i (GPU support)
      standard,           ///< Standard
      extended,           ///< Extended classical (GPU support)
      directBAMG,         ///< Direct with separation of weights (GPU support)
      modifiedExtended,   ///< Modularized extended classical (GPU support)
      modifiedExtendedI,  ///< Modularized extended+i (GPU support)
      modifiedExtendedE   ///< Modularized extended+e (GPU support)
    };

    /// AMG interpolation type for aggressive coarsening levels (HYPRE only)
    enum class AggInterpType : integer
    {
      default_,           ///< Use LAI's default option
      extendedIStage2,    ///< Extended+i 2-stage (GPU support)
      standardStage2,     ///< Standard 2-stage
      extendedStage2,     ///< Extended 2-stage (GPU support)
      multipass,          ///< Multipass (GPU support)
      modifiedExtended,   ///< Modularized Extended (GPU support)
      modifiedExtendedI,  ///< Modularized Extended+i (GPU support)
      modifiedExtendedE,  ///< Modularized Extended+e (GPU support)
      modifiedMultipass   ///< Modularized Multipass (GPU support)
    };

    /// Null space type
    enum class NullSpaceType : integer
    {
      constantModes,      ///< Constant modes
      rigidBodyModes      ///< Rigid body modes
    };

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
    CoarseningType coarseningType = CoarseningType::PMIS;           ///< Coarsening algorithm (GPUs)
    SmootherType smootherType = SmootherType::l1jacobi;             ///< Smoother type (GPUs)
#else
    CoarseningType coarseningType = CoarseningType::HMIS;           ///< Coarsening algorithm (CPUs)
    SmootherType smootherType = SmootherType::l1sgs;                ///< Smoother type (CPUs)
#endif

    integer maxLevels = 20;                                         ///< Maximum number of coarsening levels
    CycleType cycleType = CycleType::V;                             ///< AMG cycle type
    CoarseType coarseType = CoarseType::direct;                     ///< Coarse-level solver/smoother
    InterpType interpolationType = InterpType::extendedI;           ///< Interpolation algorithm
    integer interpolationMaxNonZeros = 4;                           ///< Interpolation - Max. nonzeros/row
    real64 relaxWeight = 1.0;                                       ///< Relaxation weight
    integer numSweeps = 1;                                          ///< Number of smoother sweeps
    integer numFunctions = 1;                                       ///< Number of amg functions
    integer aggressiveNumPaths = 1;                                 ///< Number of paths agg. coarsening.
    integer aggressiveNumLevels = 0;                                ///< Number of levels for aggressive coarsening.
    AggInterpType aggressiveInterpType = AggInterpType::multipass;  ///< Interp. type for agg. coarsening.
    integer aggressiveInterpMaxNonZeros = 16;                       ///< Aggressive Interpolation - Max. nonzeros/row.
    PreOrPost preOrPostSmoothing = PreOrPost::both;                 ///< Pre and/or post smoothing
    real64 threshold = 0.0;                                         ///< Threshold for "strong connections" (for classical
                                                                    ///< and smoothed-aggregation AMG)
    integer separateComponents = false;                             ///< Apply a separate component filter before AMG construction
    NullSpaceType nullSpaceType = NullSpaceType::constantModes;     ///< Null space type [constantModes,rigidBodyModes]
  }
  amg;                                                              ///< Algebraic Multigrid (AMG) parameters

  /// Multigrid reduction parameters
  struct MGR
  {
    /**
     * @brief MGR available strategies
     */
    enum class StrategyType : integer
    {
      invalid,                                     ///< default value, to ensure solver sets something
      singlePhaseReservoirFVM,                     ///< finite volume single-phase flow with wells
      singlePhaseHybridFVM,                        ///< hybrid finite volume single-phase flow
      singlePhaseReservoirHybridFVM,               ///< hybrid finite volume single-phase flow with wells
      singlePhasePoromechanics,                    ///< single phase poromechanics with finite volume single phase flow
      thermalSinglePhasePoromechanics,             ///< thermal single phase poromechanics with finite volume single phase flow
      hybridSinglePhasePoromechanics,              ///< single phase poromechanics with hybrid finite volume single phase flow
      singlePhasePoromechanicsEmbeddedFractures,   ///< single phase poromechanics with FV embedded fractures
      singlePhasePoromechanicsConformingFractures, ///< single phase poromechanics with FV conforming  fractures
      singlePhasePoromechanicsReservoirFVM,        ///< single phase poromechanics with finite volume single phase flow with wells
      compositionalMultiphaseFVM,                  ///< finite volume compositional multiphase flow
      compositionalMultiphaseHybridFVM,            ///< hybrid finite volume compositional multiphase flow
      compositionalMultiphaseReservoirFVM,         ///< finite volume compositional multiphase flow with wells
      compositionalMultiphaseReservoirHybridFVM,   ///< hybrid finite volume compositional multiphase flow with wells
      reactiveCompositionalMultiphaseOBL,          ///< finite volume reactive compositional flow with OBL
      thermalCompositionalMultiphaseFVM,           ///< finite volume thermal compositional multiphase flow
      multiphasePoromechanics,                     ///< multiphase poromechanics with finite volume compositional multiphase flow
      multiphasePoromechanicsReservoirFVM,         ///< multiphase poromechanics with finite volume compositional multiphase flow with wells
      thermalMultiphasePoromechanics,              ///< thermal multiphase poromechanics with finite volume compositional multiphase flow
      hydrofracture,                               ///< hydrofracture
      lagrangianContactMechanics,                  ///< Lagrangian contact mechanics
      solidMechanicsEmbeddedFractures              ///< Embedded fractures mechanics
    };

    StrategyType strategy = StrategyType::invalid; ///< Predefined MGR solution strategy (solver specific)
    integer separateComponents = false;            ///< Apply a separate displacement component (SDC) filter before AMG construction
    string displacementFieldName;                  ///< Displacement field name need for SDC filter
    integer areWellsShut = false;                   ///< Flag to let MGR know that wells are shut, and that jacobi can be applied to the
                                                    ///< well block
  }
  mgr;                                             ///< Multigrid reduction (MGR) parameters

  /// Incomplete factorization parameters
  struct IFact
  {
    integer fill = 0;        ///< Fill level
    real64 threshold = 0.0;  ///< Dropping threshold
  }
  ifact;                       ///< Incomplete factorization parameter struct

  /// Domain decomposition parameters
  struct DD
  {
    integer overlap = 0;   ///< Ghost overlap
  }
  dd;                      ///< Domain decomposition parameter struct
};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::SolverType,
              "direct",
              "cg",
              "gmres",
              "fgmres",
              "bicgstab",
              "preconditioner" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::PreconditionerType,
              "none",
              "jacobi",
              "l1jacobi",
              "fgs",
              "sgs",
              "l1sgs",
              "chebyshev",
              "iluk",
              "ilut",
              "icc",
              "ict",
              "amg",
              "mgr",
              "block",
              "direct",
              "bgs" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::Direct::ColPerm,
              "none",
              "MMD_AtplusA",
              "MMD_AtA",
              "colAMD",
              "metis",
              "parmetis" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::Direct::RowPerm,
              "none",
              "mc64" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::MGR::StrategyType,
              "invalid",
              "singlePhaseReservoirFVM",
              "singlePhaseHybridFVM",
              "singlePhaseReservoirHybridFVM",
              "singlePhasePoromechanics",
              "thermalSinglePhasePoromechanics",
              "hybridSinglePhasePoromechanics",
              "singlePhasePoromechanicsEmbeddedFractures",
              "singlePhasePoromechanicsConformingFractures",
              "singlePhasePoromechanicsReservoirFVM",
              "compositionalMultiphaseFVM",
              "compositionalMultiphaseHybridFVM",
              "compositionalMultiphaseReservoirFVM",
              "compositionalMultiphaseReservoirHybridFVM",
              "reactiveCompositionalMultiphaseOBL",
              "thermalCompositionalMultiphaseFVM",
              "multiphasePoromechanics",
              "multiphasePoromechanicsReservoirFVM",
              "thermalMultiphasePoromechanics",
              "hydrofracture",
              "lagrangianContactMechanics",
              "solidMechanicsEmbeddedFractures" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::CycleType,
              "V",
              "W" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::PreOrPost,
              "pre",
              "post",
              "both" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::SmootherType,
              "default",
              "jacobi",
              "l1jacobi",
              "fgs",
              "bgs",
              "sgs",
              "l1sgs",
              "chebyshev",
              "ilu0",
              "ilut",
              "ic0",
              "ict" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::CoarseType,
              "default",
              "jacobi",
              "l1jacobi",
              "fgs",
              "sgs",
              "l1sgs",
              "chebyshev",
              "direct",
              "bgs" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::CoarseningType,
              "default",
              "CLJP",
              "RugeStueben",
              "Falgout",
              "PMIS",
              "HMIS" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::InterpType,
              "default",
              "modifiedClassical",
              "direct",
              "multipass",
              "extendedI",
              "standard",
              "extended",
              "directBAMG",
              "modifiedExtended",
              "modifiedExtendedI",
              "modifiedExtendedE" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::AggInterpType,
              "default",
              "extendedIStage2",
              "standardStage2",
              "extendedStage2",
              "multipass",
              "modifiedExtended",
              "modifiedExtendedI",
              "modifiedExtendedE",
              "modifiedMultipass" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::NullSpaceType,
              "constantModes",
              "rigidBodyModes" );

} /* namespace geos */

#endif /*GEOS_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_ */
