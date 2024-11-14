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
 * @file HypreMGR.cpp
 */

#include "HypreMGR.hpp"


#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseReservoirHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/HybridSinglePhasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/Hydrofracture.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ImmiscibleMultiphaseFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/LagrangianContactMechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/MultiphasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/MultiphasePoromechanicsReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ReactiveCompositionalMultiphaseOBL.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhaseHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsEmbeddedFractures.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsConformingFractures.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhaseReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhaseReservoirHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalCompositionalMultiphaseFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalCompositionalMultiphaseReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalSinglePhasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalMultiphasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SolidMechanicsEmbeddedFractures.hpp"

#include "LvArray/src/output.hpp"

namespace geos
{

void hypre::mgr::createMGR( LinearSolverParameters const & params,
                            DofManager const * const dofManager,
                            HyprePrecWrapper & precond,
                            HypreMGRData & mgrData )
{
  GEOS_ERROR_IF( dofManager == nullptr, "MGR preconditioner requires a DofManager instance" );

  GEOS_LAI_CHECK_ERROR( HYPRE_MGRCreate( &precond.ptr ) );

  // Hypre's parameters to use MGR as a preconditioner
  GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetTol( precond.ptr, 0.0 ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetMaxIter( precond.ptr, 1 ) );

  // Disabling Option 2 (0x2): In this context, we use MGR as a preconditioner for a Krylov method (e.g., GMRES),
  // and our interest lies in the Krylov method's convergence history, not MGR's. Hence, we turn off the second bit
  // (0x2) of logLevel. For detailed logLevel codes, see HYPRE_MGRSetPrintLevel documentation. Additionally,
  // we subtract one from the input logLevel value because "1" is reserved in GEOS for setup and solve time logging.
  HYPRE_Int logLevel = LvArray::math::max( LvArray::integerConversion< HYPRE_Int >( params.logLevel - 1 ), LvArray::integerConversion< HYPRE_Int >( 0 ) );
  logLevel &= ~0x2;
  GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetPrintLevel( precond.ptr, logLevel ) );

  array1d< int > const numComponentsPerField = dofManager->numComponentsPerField();
  dofManager->getLocalDofComponentLabels( mgrData.pointMarkers );

  if( params.logLevel >= 1 )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "        MGR preconditioner: numComponentsPerField = {}", numComponentsPerField ) );
  }
  if( params.logLevel >= 1024 )
  {
    GEOS_LOG_RANK( GEOS_FMT( "        MGR preconditioner: pointMarkers = {}", mgrData.pointMarkers ) );
  }

  switch( params.mgr.strategy )
  {
    case LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseFVM:
    {
      setStrategy< CompositionalMultiphaseFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseHybridFVM:
    {
      setStrategy< CompositionalMultiphaseHybridFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseReservoirFVM:
    {
      setStrategy< CompositionalMultiphaseReservoirFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseReservoirHybridFVM:
    {
      setStrategy< CompositionalMultiphaseReservoirHybridFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::reactiveCompositionalMultiphaseOBL:
    {
      setStrategy< ReactiveCompositionalMultiphaseOBL >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::thermalCompositionalMultiphaseFVM:
    {
      setStrategy< ThermalCompositionalMultiphaseFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::thermalCompositionalMultiphaseReservoirFVM:
    {
      setStrategy< ThermalCompositionalMultiphaseReservoirFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::hybridSinglePhasePoromechanics:
    {
      setStrategy< HybridSinglePhasePoromechanics >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::thermalSinglePhasePoromechanics:
    {
      setStrategy< ThermalSinglePhasePoromechanics >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::hydrofracture:
    {
      setStrategy< Hydrofracture >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::lagrangianContactMechanics:
    {
      setStrategy< LagrangianContactMechanics >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::multiphasePoromechanics:
    {
      setStrategy< MultiphasePoromechanics >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::multiphasePoromechanicsReservoirFVM:
    {
      setStrategy< MultiphasePoromechanicsReservoirFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::thermalMultiphasePoromechanics:
    {
      setStrategy< ThermalMultiphasePoromechanics >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM:
    {
      setStrategy< SinglePhaseHybridFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics:
    {
      setStrategy< SinglePhasePoromechanics >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsEmbeddedFractures:
    {
      setStrategy< SinglePhasePoromechanicsEmbeddedFractures >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsConformingFractures:
    {
      setStrategy< SinglePhasePoromechanicsConformingFractures >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsReservoirFVM:
    {
      setStrategy< SinglePhasePoromechanicsReservoirFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::singlePhaseReservoirFVM:
    {
      setStrategy< SinglePhaseReservoirFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::singlePhaseReservoirHybridFVM:
    {
      setStrategy< SinglePhaseReservoirHybridFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::solidMechanicsEmbeddedFractures:
    {
      setStrategy< SolidMechanicsEmbeddedFractures >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::immiscibleMultiphaseFVM:
    {
      setStrategy< ImmiscibleMultiphaseFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    default:
    {
      GEOS_ERROR( "Unsupported MGR strategy: " << params.mgr.strategy );
    }
  }

  GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseSolver( precond.ptr,
                                                  mgrData.coarseSolver.solve,
                                                  mgrData.coarseSolver.setup,
                                                  mgrData.coarseSolver.ptr ) );
  precond.setup = HYPRE_MGRSetup;
  precond.solve = HYPRE_MGRSolve;
  precond.destroy = HYPRE_MGRDestroy;
}

} // namespace geos
