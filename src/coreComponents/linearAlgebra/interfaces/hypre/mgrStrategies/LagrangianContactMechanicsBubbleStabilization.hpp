/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LagrangianContactMechanicsBubbleStabilization.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRLAGRANGIACONTACTMECHANICSBUBBLESTAB_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRLAGRANGIACONTACTMECHANICSBUBBLESTAB_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief LagrangianContactMechanicsBubbleStabilization strategy
 *
 * Contact mechanics with face-centered lagrangian multipliers
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = displacement bubble function, x-component
 * dofLabel: 4 = displacement bubble function, y-component
 * dofLabel: 5 = displacement bubble function, z-component
 * dofLabel: 6 = face-centered lagrange multiplier (tn)
 * dofLabel: 7 = face-centered lagrange multiplier (tt1)
 * dofLabel: 8 = face-centered lagrange multiplier (tt2)
 *
 * 2-level MGR strategy:
 * Level 0:
 * 1. F-points: bubbles (3,4,5), C-points: displacements + multipliers (0,1,2,6,7,8)
 * 2. F-points smoother: l1jacobi
 * 3. Global smoother: none
 * Level 1:
 * 1. F-points: multipliers (6,7,8), C-points: displacements (0,1,2)
 * 2. F-points smoother: l1jacobi
 * 3. Global smoother: none
 *
 * C-points coarse-grid/Schur complement solver: boomer AMG

 */
class LagrangianContactMechanicsBubbleStabilization : public MGRStrategyBase< 2 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit LagrangianContactMechanicsBubbleStabilization( arrayView1d< int const > const & )
    : MGRStrategyBase( 9 )
  {
    // Level 0: we keep all three displacements and the Lagrange Multipliers
    m_labels[0] = { 0, 1, 2, 6, 7, 8 };
    // Level 1: we keep all three displacements
    m_labels[1] = { 0, 1, 2 };

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]          = MGRFRelaxationType::l1jacobi;
    m_levelFRelaxIters[0]         = 1;
    m_levelInterpType[0]          = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;

    // Level 1
    m_levelFRelaxType[1]          = MGRFRelaxationType::l1jacobi;
    m_levelFRelaxIters[1]         = 1;
    m_levelInterpType[1]          = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[1]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1]  = MGRGlobalSmootherType::none;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param mgrParams MGR configuration parameters
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const & mgrParams,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    setReduction( precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the displacement reduced system
    setDisplacementAMG( mgrData.coarseSolver, mgrParams.separateComponents );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRLAGRANGIACONTACTMECHANICSBUBBLESTAB_HPP_*/
