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
 * @file HyprePreconditioner.cpp
 */

#include "HyprePreconditioner.hpp"
#include "HypreMGR.hpp"

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <_hypre_utilities.h>
#include <_hypre_parcsr_ls.h>

#include <_hypre_IJ_mv.h>

#include <cfenv>

namespace geos
{

/**
 * @brief Container for hypre preconditioner null space data
 */
struct HypreNullSpace
{
  array1d< HYPRE_ParVector > vectors; ///< Hypre vectors containing the near null kernel
};

namespace
{

void convertRigidBodyModes( arrayView1d< HypreVector > const & nearNullKernel,
                            array1d< HYPRE_ParVector > & nullSpacePointer )
{
  if( nearNullKernel.empty() )
  {
    return;
  }
  else
  {
    localIndex dim = 0;
    if( nearNullKernel.size() == 3 )
    {
      dim = 2;
    }
    else if( nearNullKernel.size() == 6 )
    {
      dim = 3;
    }
    else
    {
      GEOS_ERROR( "Hypre preconditioner: rigid body modes can be either 3 or 6. Current number: " << nearNullKernel.size() );
    }
    localIndex const numRotations = LvArray::integerConversion< HYPRE_Int >( nearNullKernel.size() - dim );
    nullSpacePointer.resize( numRotations );
    for( localIndex k = 0; k < numRotations; ++k )
    {
      nullSpacePointer[k] = nearNullKernel[dim+k].unwrapped();
    }
  }
}

void createAMG( LinearSolverParameters const & params,
                HypreNullSpace const & nullSpace,
                HyprePrecWrapper & precond )
{
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &precond.ptr ) );

  // Hypre's parameters to use BoomerAMG as a preconditioner
  HYPRE_Int logLevel = (params.logLevel == 2 || params.logLevel >= 4) ? 1 : 0;

  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( precond.ptr, 0.0 ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( precond.ptr, 1 ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( precond.ptr, logLevel ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( precond.ptr, params.dofsPerNode ) );

  // Set maximum number of multigrid levels (default 25)
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxLevels( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.amg.maxLevels ) ) );

  // Set type of cycle (1: V-cycle (default); 2: W-cycle)
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleType( precond.ptr, hypre::getAMGCycleType( params.amg.cycleType ) ) );

  if( params.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes && !nullSpace.vectors.empty() )
  {
    // Set of options used in MFEM
    // Nodal coarsening options (nodal coarsening is required for this solver)
    // See hypre's new_ij driver and the paper for descriptions.

    // For further information, see:
    // Improving algebraic multigrid interpolation operators for linear elasticity problems
    // A. H. Baker Tz. V. Kolev U. M. Yang
    // Numerical Linear Algebra with Applications (2010) 17 (2-3), 495-517
    // doi:10.1002/nla.688

    HYPRE_Int const nodal                 = 4; // strength reduction norm: 1, 3 or 4
    HYPRE_Int const nodal_diag            = 1; // diagonal in strength matrix: 0, 1 or 2
    HYPRE_Int const relax_coarse          = 8; // smoother on the coarsest grid: 8, 99 or 29

    // Elasticity interpolation options
    HYPRE_Int const interp_vec_variant    = 2; // 1 = GM-1, 2 = GM-2, 3 = LN
    HYPRE_Int const q_max                 = 4; // max elements per row for each Q
    HYPRE_Int const smooth_interp_vectors = 1; // smooth the rigid-body modes?

    // Optionally pre-process the interpolation matrix through iterative weight
    // refinement (this is generally applicable for any system)
    HYPRE_Int const interp_refine         = 1;

    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNodal( precond.ptr, nodal ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNodalDiag( precond.ptr, nodal_diag ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleRelaxType( precond.ptr, relax_coarse, 3 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVecVariant( precond.ptr, interp_vec_variant ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVecQMax( precond.ptr, q_max ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothInterpVectors( precond.ptr, smooth_interp_vectors ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpRefine( precond.ptr, interp_refine ) );

    // Add user-defined null space / rigid body mode support
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVectors( precond.ptr, nullSpace.vectors.size(), nullSpace.vectors.data() ) );
  }

  // Set smoother to be used (other options available, see hypre's documentation)
  // (default "gaussSeidel", i.e. local symmetric Gauss-Seidel)

  if( params.amg.smootherType == LinearSolverParameters::AMG::SmootherType::ilu0 ||
      params.amg.smootherType == LinearSolverParameters::AMG::SmootherType::ilut )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothType( precond.ptr, 5 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_ILUSetType( precond.ptr, hypre::getILUType( params.amg.smootherType ) ) );
  }
  else
  {
    HYPRE_Int const relaxType = hypre::getAMGRelaxationType( params.amg.smootherType );
    if( relaxType >= 0 )
    {
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( precond.ptr, relaxType ) );
    }

    if( params.amg.smootherType == LinearSolverParameters::AMG::SmootherType::chebyshev )
    {
      // Set order for Chebyshev smoother valid options 1, 2 (default), 3, 4)
      if( ( 0 < params.amg.numSweeps ) && ( params.amg.numSweeps < 5 ) )
      {
        GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetChebyOrder( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.amg.numSweeps ) ) );
      }
    }
  }

  // Set relaxation weight
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxWt( precond.ptr, params.amg.relaxWeight ) );

  // Coarsening options: Only PMIS is supported on GPU
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( precond.ptr, hypre::getAMGCoarseningType( params.amg.coarseningType ) ) );

  // Interpolation options
  {
    HYPRE_Int const interpType = hypre::getAMGInterpolationType( params.amg.interpolationType );

    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpType( precond.ptr, interpType ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPMaxElmts( precond.ptr, params.amg.interpolationMaxNonZeros ) );
  }

  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( precond.ptr, params.amg.numFunctions ) );

  if( params.amg.aggressiveNumLevels )
  {
    HYPRE_Int const aggInterpType = hypre::getAMGAggressiveInterpolationType( params.amg.aggressiveInterpType );

    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumPaths( precond.ptr, params.amg.aggressiveNumPaths ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggNumLevels( precond.ptr, params.amg.aggressiveNumLevels ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggInterpType( precond.ptr, aggInterpType ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggPMaxElmts( precond.ptr, params.amg.aggressiveInterpMaxNonZeros ) );
  }

  // Set coarsest level solver
  HYPRE_Int const coarseType = hypre::getAMGCoarseType( params.amg.coarseType );
  if( coarseType >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleRelaxType( precond.ptr, coarseType, 3 ) );
  }

  // Set the number of sweeps
  switch( params.amg.preOrPostSmoothing )
  {
    case LinearSolverParameters::AMG::PreOrPost::both:
    {
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.amg.numSweeps ) ) );
      break;
    }
    case LinearSolverParameters::AMG::PreOrPost::pre:
    {
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.amg.numSweeps ), 1 ) );
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( precond.ptr, 0, 2 ) );
      break;
    }
    case LinearSolverParameters::AMG::PreOrPost::post:
    {
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( precond.ptr, 0, 1 ) );
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.amg.numSweeps ), 2 ) );
      break;
    }
    default:
    {
      GEOS_ERROR( "Unsupported Hypre/AMG smoothing direction: " << params.amg.preOrPostSmoothing );
    }
  }

  // Set strength of connection
  if( params.amg.threshold > 0.0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( precond.ptr, params.amg.threshold ) );
  }

  precond.setup = HYPRE_BoomerAMGSetup;
  precond.solve = HYPRE_BoomerAMGSolve;
  precond.destroy = HYPRE_BoomerAMGDestroy;
}

void createILU( LinearSolverParameters const & params,
                HyprePrecWrapper & precond )
{
  GEOS_LAI_CHECK_ERROR( HYPRE_ILUCreate( &precond.ptr ) );

  // Hypre's parameters to use ParCSR ILU as a preconditioner
  GEOS_LAI_CHECK_ERROR( HYPRE_ILUSetMaxIter( precond.ptr, 1 ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_ILUSetTol( precond.ptr, 0.0 ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_ILUSetType( precond.ptr, hypre::getILUType( params.preconditionerType ) ) );

  if( params.ifact.fill >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.ifact.fill ) ) );
  }
  if( params.ifact.threshold >= 0 && params.preconditionerType == LinearSolverParameters::PreconditionerType::ilut )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_ILUSetDropThreshold( precond.ptr, params.ifact.threshold ) );
  }

  precond.setup = HYPRE_ILUSetup;
  precond.solve = HYPRE_ILUSolve;
  precond.destroy = HYPRE_ILUDestroy;
}

void createRelaxation( LinearSolverParameters const & params,
                       HyprePrecWrapper & precond )
{
  GEOS_LAI_CHECK_ERROR( hypre::relaxationCreate( precond.ptr, hypre::getRelaxationType( params.preconditionerType ) ) );
  precond.setup = hypre::relaxationSetup;
  precond.solve = hypre::relaxationSolve;
  precond.destroy = hypre::relaxationDestroy;
}

} // namespace

HyprePreconditioner::HyprePreconditioner( LinearSolverParameters params )
  : Base{},
  m_params( std::move( params ) ),
  m_nullSpace( std::make_unique< HypreNullSpace >() )
{}

HyprePreconditioner::HyprePreconditioner( LinearSolverParameters params,
                                          arrayView1d< HypreVector > const & nearNullKernel )
  : HyprePreconditioner( std::move( params ) )
{
  if( m_params.preconditionerType == LinearSolverParameters::PreconditionerType::amg &&
      m_params.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes )
  {
    convertRigidBodyModes( nearNullKernel, m_nullSpace->vectors );
  }
}

HyprePreconditioner::~HyprePreconditioner()
{
  HyprePreconditioner::clear();
}

void HyprePreconditioner::create( DofManager const * const dofManager )
{
  switch( m_params.preconditionerType )
  {
    case LinearSolverParameters::PreconditionerType::none:
    {
      m_precond->setup = (HyprePrecWrapper::SetupFunc) hypre_ParKrylovIdentitySetup;
      m_precond->solve = (HyprePrecWrapper::SetupFunc) hypre_ParKrylovIdentity;
      break;
    }
    case LinearSolverParameters::PreconditionerType::jacobi:
    case LinearSolverParameters::PreconditionerType::fgs:
    case LinearSolverParameters::PreconditionerType::bgs:
    case LinearSolverParameters::PreconditionerType::sgs:
    case LinearSolverParameters::PreconditionerType::l1jacobi:
    case LinearSolverParameters::PreconditionerType::chebyshev:
    case LinearSolverParameters::PreconditionerType::l1sgs:
    {
      createRelaxation( m_params, *m_precond );
      break;
    }
    case LinearSolverParameters::PreconditionerType::amg:
    {
      createAMG( m_params, *m_nullSpace, *m_precond );
      break;
    }
    case LinearSolverParameters::PreconditionerType::mgr:
    {
      m_mgrData = std::make_unique< HypreMGRData >();
      hypre::mgr::createMGR( m_params, dofManager, *m_precond, *m_mgrData );
      break;
    }
    case LinearSolverParameters::PreconditionerType::iluk:
    case LinearSolverParameters::PreconditionerType::ilut:
    {
      createILU( m_params, *m_precond );
      break;
    }
    case LinearSolverParameters::PreconditionerType::direct:
    {
      m_precond->solve = hypre::SuperLUDistSolve;
      m_precond->destroy = hypre::SuperLUDistDestroy;
      break;
    }
    default:
    {
      GEOS_ERROR( "Preconditioner type not supported in hypre interface: " << m_params.preconditionerType );
    }
  }
}

HypreMatrix const & HyprePreconditioner::setupPreconditioningMatrix( HypreMatrix const & mat )
{
  GEOS_MARK_FUNCTION;

  if( m_params.preconditionerType == LinearSolverParameters::PreconditionerType::mgr && m_params.mgr.separateComponents )
  {
    GEOS_LAI_ASSERT_MSG( mat.dofManager() != nullptr, "MGR preconditioner requires a DofManager instance" );
    HypreMatrix Pu;
    HypreMatrix Auu;
    {
      Stopwatch timer( m_makeRestrictorTime );
      mat.dofManager()->makeRestrictor( { { m_params.mgr.displacementFieldName, { 3, true } } }, mat.comm(), true, Pu );
    }
    {
      Stopwatch timer( m_computeAuuTime );
      mat.multiplyPtAP( Pu, Auu );
    }
    {
      Stopwatch timer( m_componentFilterTime );
      Auu.separateComponentFilter( m_precondMatrix, m_params.dofsPerNode );
    }
  }
  else if( m_params.preconditionerType == LinearSolverParameters::PreconditionerType::amg && m_params.amg.separateComponents )
  {
    Stopwatch timer( m_componentFilterTime );
    mat.separateComponentFilter( m_precondMatrix, m_params.dofsPerNode );
    return m_precondMatrix;
  }
  return mat;
}

void HyprePreconditioner::setup( Matrix const & mat )
{
  GEOS_MARK_FUNCTION;

  if( !m_precond )
  {
    m_precond = std::make_unique< HyprePrecWrapper >();
    create( mat.dofManager() );
  }
  HypreMatrix const & precondMat = setupPreconditioningMatrix( mat );
  Base::setup( precondMat );

  // To be able to use Hypre preconditioner (e.g., BoomerAMG) we need to disable floating point exceptions
  {
    LvArray::system::FloatingPointExceptionGuard guard( FE_ALL_EXCEPT );

    // Perform setup of the MGR mechanics F-solver with SDC matrix, if used
    if( m_mgrData && m_mgrData->mechSolver.ptr && m_mgrData->mechSolver.setup )
    {
//      GEOS_LAI_CHECK_ERROR( m_mgrData->mechSolver.setup( m_mgrData->mechSolver.ptr, m_precondMatrix.unwrapped(), nullptr, nullptr ) );
    }

    // Perform setup of the main solver, if needed
    if( m_precond->setup )
    {
      HYPRE_Int const ierr = m_precond->setup( m_precond->ptr, precondMat.unwrapped(), nullptr, nullptr );
      if( ierr != 0 )
      {
        if( m_params.logLevel > 0 )
        {
          HYPRE_PrintErrorMessages( precondMat.comm() );
          GEOS_ERROR( "HyprePreconditioner setup failed" );
        }
        else
        {
          GEOS_ERROR( "HyprePreconditioner setup failed, use logLevel > 0 to get more information" );
        }
      }
    }
    else if( m_params.preconditionerType == LinearSolverParameters::PreconditionerType::direct )
    {
      // Special handling for hypre's SuperLU_Dist interface: it combines Create and Setup methods in one,
      // and thus we have to reallocate the entire solver data structure
      if( m_precond->ptr && m_precond->destroy )
      {
        m_precond->destroy( m_precond->ptr );
      }
#if defined(GEOS_USE_SUPERLU_DIST)
      hypre_SLUDistSetup( &m_precond->ptr, precondMat.unwrapped(), 0 );
#endif
    }
  }
}

void HyprePreconditioner::apply( Vector const & src,
                                 Vector & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT( dst.ready() );
  GEOS_LAI_ASSERT_EQ( src.localSize(), numLocalCols() );
  GEOS_LAI_ASSERT_EQ( dst.localSize(), numLocalRows() );

  // Needed to avoid accumulation inside HYPRE solver phase
  dst.zero();

  GEOS_LAI_CHECK_ERROR( m_precond->solve( m_precond->ptr, matrix().unwrapped(), src.unwrapped(), dst.unwrapped() ) );
  dst.touch();
}

void HyprePreconditioner::clear()
{
  Base::clear();
  if( m_precond && m_precond->ptr && m_precond->destroy )
  {
    GEOS_LAI_CHECK_ERROR( m_precond->destroy( m_precond->ptr ) );
  }
  if( m_mgrData && m_mgrData->coarseSolver.ptr && m_mgrData->coarseSolver.destroy )
  {
    GEOS_LAI_CHECK_ERROR( m_mgrData->coarseSolver.destroy( m_mgrData->coarseSolver.ptr ) );
  }
  if( m_mgrData && m_mgrData->mechSolver.ptr && m_mgrData->mechSolver.destroy )
  {
    GEOS_LAI_CHECK_ERROR( m_mgrData->mechSolver.destroy( m_mgrData->mechSolver.ptr ) );
  }
  m_precond.reset();
  m_mgrData.reset();
}

HyprePrecWrapper const & HyprePreconditioner::unwrapped() const
{
  GEOS_LAI_ASSERT( m_precond );
  return *m_precond;
}

}
