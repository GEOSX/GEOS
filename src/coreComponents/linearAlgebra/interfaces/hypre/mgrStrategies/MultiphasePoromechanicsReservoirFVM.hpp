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
 * @file MultiphasePoromechanicsReservoirFVM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEPOROMECHANICSRESERVOIRFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEPOROMECHANICSRESERVOIRFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief MultiphasePoromechanicsReservoirFVM strategy.
 *
 * Labels description stored in point_marker_array
 *   - dofLabel: 0                = nodal displacement, x-component
 *   - dofLabel: 1                = nodal displacement, y-component
 *   - dofLabel: 2                = nodal displacement, z-component
 *   - dofLabel: 3                = pressure
 *   - dofLabel: 4                = density
 *             ...                = densities
 *   - dofLabel: numResLabels + 2 = density
 *   - dofLabel: numResLabels + 3 = well pressure
 *   - dofLabel: numResLabels + 4 = well density
 *             ...                = well densities
 *   - dofLabel: numResLabels + numWellLabels + 1 = well density
 *   - dofLabel: numResLabels + numWellLabels + 2 = well rate
 *
 * 4-level MGR reduction strategy
 *   - 1st level: eliminate displacements (0,1,2)
 *   - 2nd level: eliminate the well block
 *   - 3nd level: eliminate the reservoir density associated with the volume constraint
 *   - 4nd level: eliminate the other reservoir densities
 *   - The coarse grid is solved with BoomerAMG.
 *
 */
class MultiphasePoromechanicsReservoirFVM : public MGRStrategyBase< 4 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit MultiphasePoromechanicsReservoirFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] + numComponentsPerField[2] ) )
  {
    // Level 0: eliminate displacement degrees of freedom
    m_labels[0].resize( m_numBlocks - 3 );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 3 );

    HYPRE_Int const numResLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[1] );

    // Level 1: eliminate the well block
    m_labels[1].resize( numResLabels );
    std::iota( m_labels[1].begin(), m_labels[1].end(), 3 );
    // Level 2: eliminate reservoir last density which corresponds to the volume constraint equation
    m_labels[2].resize( numResLabels - 1 );
    std::iota( m_labels[2].begin(), m_labels[2].end(), 3 );
    // Level 3: eliminate the remaining reservoir densities
    m_labels[3].push_back( 3 );

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]          = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[0]         = 1;
    m_levelInterpType[0]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::nonGalerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;

    // Level 1
    m_levelFRelaxType[1]          = MGRFRelaxationType::gsElimWInverse;
    m_levelFRelaxIters[1]         = 1;
    m_levelInterpType[1]          = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[1]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1]  = MGRGlobalSmootherType::none;

    // Level 2
    m_levelFRelaxType[2]          = MGRFRelaxationType::jacobi; //default, i.e. Jacobi
    m_levelFRelaxIters[2]         = 1;
    m_levelInterpType[2]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[2]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[2]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[2]  = MGRGlobalSmootherType::none;

    // Level 3
    m_levelFRelaxType[3]          = MGRFRelaxationType::none;
    m_levelInterpType[3]          = MGRInterpolationType::injection;
    m_levelRestrictType[3]        = MGRRestrictionType::blockColLumped; // True-IMPES
    m_levelCoarseGridMethod[3]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[3]  = MGRGlobalSmootherType::ilu0;
    m_levelGlobalSmootherIters[3] = 1;
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
      m_levelFRelaxType[1] = MGRFRelaxationType::jacobi;
    }

    setReduction( precond, mgrData );

    // Configure the BoomerAMG solver used as F-relaxation for the first level
    setMechanicsFSolver( precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEPOROMECHANICSRESERVOIRFVM_HPP_*/
