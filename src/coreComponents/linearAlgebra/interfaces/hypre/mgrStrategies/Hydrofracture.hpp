/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Hydrofracture.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRHYDROFRACTURE_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRHYDROFRACTURE_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief Hydrofracture strategy.
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = pressure
 *
 * Ingredients:
 * 1. F-points displacement (0,1,2), C-points pressure (3)
 * 2. F-points smoother: AMG, single V-cycle, separate displacement components
 * 3. C-points coarse-grid/Schur complement solver: boomer AMG
 * 4. Global smoother: none
 */
class Hydrofracture : public MGRStrategyBase< 1 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit Hydrofracture( arrayView1d< int const > const & )
    : MGRStrategyBase( 4 )
  {
    m_labels[0].push_back( 3 );

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]          = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[0]         = 1;
    m_levelInterpType[0]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;
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

    // Configure the BoomerAMG solver used as F-relaxation for the first level
    setMechanicsFSolver( precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMinCoarseSize( mgrData.coarseSolver.ptr, 1000 ) );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRHYDROFRACTURE_HPP_*/
