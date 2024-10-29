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
 * @file SolidMechanicsEmbeddedFractures.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSOLIDMECHANICSEMBEDDEDFRACTURES_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSOLIDMECHANICSEMBEDDEDFRACTURES_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief SolidmechanicsEmbeddedFractures strategy.
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = wn
 * dofLabel: 4 = wt1
 * dofLabel: 5 = wt2
 *
 * Ingredients:
 * 1. F-points w (3 - 4 - 5), C-points displacement (0 - 1 - 2)
 * 2. F-points smoother: gaussian elimination
 * 3. C-points coarse-grid/Schur complement solver: boomer AMG
 * 4. Global smoother: none
 */
class SolidMechanicsEmbeddedFractures : public MGRStrategyBase< 1 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit SolidMechanicsEmbeddedFractures( arrayView1d< int const > const & )
    : MGRStrategyBase( 6 )
  {
    // we keep u
    m_labels[0].push_back( 0 );
    m_labels[0].push_back( 1 );
    m_labels[0].push_back( 2 );

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]         = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[0]        = 1;
    m_levelInterpType[0]         = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0] = MGRGlobalSmootherType::none;
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

    // Configure the BoomerAMG solver used as mgr coarse solver for the displacement reduced system
    setDisplacementAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSOLIDMECHANICSEMBEDDEDFRACTURES_HPP_*/
