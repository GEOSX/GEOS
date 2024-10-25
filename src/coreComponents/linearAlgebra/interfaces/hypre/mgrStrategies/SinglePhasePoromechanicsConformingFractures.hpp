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
 * @file SinglePhasePoromechanicsConformingFractures.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief SinglePhasePoromechanicsConformingFractures strategy.
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = face-centered lagrange multiplier (tn)
 * dofLabel: 4 = face-centered lagrange multiplier (tt1)
 * dofLabel: 5 = face-centered lagrange multiplier (tt2)
 * dofLabel: 6 = pressure (cell elem + fracture elems)

 *
 * Ingredients:
 * 1. Level 1: F-points displacement (3,4,5), C-points pressure (0,1,2,6)
 * 2. Level 2: F-points displacement (0,1,2), C-points pressure (6)
 * 2. F-points smoother: BoomerAMG, single V-cycle
 * 3. C-points coarse-grid/Schur complement solver: BoomerAMG
 * 4. Global smoother: none
 */
class SinglePhasePoromechanicsConformingFractures : public MGRStrategyBase< 2 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit SinglePhasePoromechanicsConformingFractures( arrayView1d< int const > const & )
    : MGRStrategyBase( 7 )
  {

    // we keep u and p
    m_labels[0].push_back( 0 );
    m_labels[0].push_back( 1 );
    m_labels[0].push_back( 2 );
    m_labels[0].push_back( 6 );
    // we keep p
    m_labels[1].push_back( 6 );

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]          = MGRFRelaxationType::none;
    m_levelFRelaxIters[0]         = 0;

    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::ilu0;
    m_levelGlobalSmootherIters[0] = 1;

    m_levelInterpType[0]          = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;

    // Level 1
    m_levelFRelaxType[1]         = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[1]        = 1;
    m_levelGlobalSmootherType[1] = MGRGlobalSmootherType::none;
    m_levelInterpType[1]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]   = MGRCoarseGridMethod::nonGalerkin;

  }

  /**
   * @brief Setup the MGR strategy.
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const &,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    setReduction( precond, mgrData );

    // Configure the BoomerAMG solver used as F-relaxation for the second level
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.mechSolver.ptr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.mechSolver.ptr, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.mechSolver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( mgrData.mechSolver.ptr, 1.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( mgrData.mechSolver.ptr, 0.6 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.mechSolver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( mgrData.mechSolver.ptr, 3 ) );

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( mgrData.mechSolver.ptr, hypre::getAMGCoarseningType( LinearSolverParameters::AMG::CoarseningType::PMIS ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.mechSolver.ptr, hypre::getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::chebyshev ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.mechSolver.ptr, 1 ) );
#else
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.mechSolver.ptr, 1 ) );
#endif
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetFSolverAtLevel( precond.ptr, mgrData.mechSolver.ptr, 1 ) );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_*/
