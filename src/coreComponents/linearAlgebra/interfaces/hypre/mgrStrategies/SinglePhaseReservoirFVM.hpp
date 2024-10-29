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
 * @file SinglePhaseReservoirFVM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief SinglePhaseReservoirFVM strategy.
 *
 * Labels description stored in point_marker_array
 *   dofLabel: 0 = reservoir pressure (numResLabels = 1)
 *   dofLabel: 1 = well pressure
 *   dofLabel: 2 = well rate (numWellLabels = 2)
 *
 * Ingredients
 *
 * 1. F-points well vars, C-points cell-centered pressure
 * 2. F-points smoother: jacobi (soon, direct solver)
 * 3. C-points coarse-grid/Schur complement solver: BoomerAMG
 */
class SinglePhaseReservoirFVM : public MGRStrategyBase< 1 >
{
public:
  /**
   * @brief Constructor.
   */
  explicit SinglePhaseReservoirFVM( arrayView1d< int const > const & )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( 3 ) )
  {
    // Level 0: eliminate the well variables, and just keep the cell-centered pressures
    m_labels[0].push_back( 0 );

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]          = MGRFRelaxationType::gsElimWInverse;
    m_levelFRelaxIters[0]        = 1;
    m_levelInterpType[0]         = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0] = MGRGlobalSmootherType::none;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param mgrParams parameters for the configuration of the MGR recipe
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const & mgrParams,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    // if the wells are shut, using Gaussian elimination as F-relaxation for the well block is an overkill
    // in that case, we just use Jacobi
    if( mgrParams.areWellsShut )
    {
      m_levelFRelaxType[0] = MGRFRelaxationType::jacobi;
    }

    setReduction( precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRFVM_HPP_*/
