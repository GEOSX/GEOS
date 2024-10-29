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
 * @file testExternalSolvers.cpp
 */

#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <gtest/gtest.h>

using namespace geos;

///////////////////////////////////////////////////////////////////////////////////////

LinearSolverParameters params_DirectSerial()
{
  LinearSolverParameters parameters;
  parameters.solverType = geos::LinearSolverParameters::SolverType::direct;
  parameters.direct.parallel = 0;
  return parameters;
}

LinearSolverParameters params_DirectParallel()
{
  LinearSolverParameters parameters;
  parameters.solverType = geos::LinearSolverParameters::SolverType::direct;
  parameters.direct.parallel = 1;
  return parameters;
}

LinearSolverParameters params_GMRES_ILU()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.solverType = LinearSolverParameters::SolverType::gmres;
  parameters.preconditionerType = LinearSolverParameters::PreconditionerType::iluk;
  parameters.ifact.fill = 1;
  return parameters;
}

LinearSolverParameters params_CG_SGS()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.isSymmetric = true;
  parameters.solverType = LinearSolverParameters::SolverType::cg;
  parameters.preconditionerType = LinearSolverParameters::PreconditionerType::sgs;
  return parameters;
}

LinearSolverParameters params_GMRES_AMG()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.solverType = LinearSolverParameters::SolverType::gmres;
  parameters.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  parameters.amg.smootherType = geos::LinearSolverParameters::AMG::SmootherType::fgs;
  parameters.amg.coarseType = geos::LinearSolverParameters::AMG::CoarseType::direct;
  return parameters;
}

LinearSolverParameters params_CG_AMG()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.isSymmetric = true;
  parameters.solverType = LinearSolverParameters::SolverType::cg;
  parameters.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  parameters.amg.smootherType = geos::LinearSolverParameters::AMG::SmootherType::sgs;
  parameters.amg.coarseType = geos::LinearSolverParameters::AMG::CoarseType::direct;
  return parameters;
}

///////////////////////////////////////////////////////////////////////////////////////

template< typename LAI >
class SolverTestBase : public ::testing::Test
{
public:

  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;

protected:

  Matrix matrix;
  real64 cond_est = 1.0;

  void test( LinearSolverParameters const & params )
  {
    // Create a random "true" solution vector
    Vector sol_true;
    sol_true.create( matrix.numLocalCols(), matrix.comm() );
    sol_true.rand( 1984 );

    // Create and compute the right-hand side vector
    Vector rhs;
    rhs.create( matrix.numLocalRows(), matrix.comm() );
    matrix.apply( sol_true, rhs );

    // Create and zero out the computed solution vector
    Vector sol_comp;
    sol_comp.create( sol_true.localSize(), sol_true.comm() );
    sol_comp.zero();

    // Create the solver and solve the system
    auto solver = LAI::createSolver( params );
    solver->setup( matrix );
    solver->solve( rhs, sol_comp );
    EXPECT_TRUE( solver->result().success() );

    // Check that solution is within epsilon of true
    Vector sol_diff( sol_comp );
    sol_diff.axpy( -1.0, sol_true );
    real64 const relTol = cond_est * params.krylov.relTolerance;
    EXPECT_LT( sol_diff.norm2() / sol_true.norm2(), relTol );
  }
};

///////////////////////////////////////////////////////////////////////////////////////

template< typename LAI >
class SolverTestLaplace2D : public SolverTestBase< LAI >
{
public:

  using Base = SolverTestBase< LAI >;
  using Matrix = typename Base::Matrix;
  using Vector = typename Base::Vector;

protected:

  void SetUp() override
  {
    globalIndex constexpr n = 100;
    geos::testing::compute2DLaplaceOperator( MPI_COMM_GEOSX, n, this->matrix );

    // Condition number for the Laplacian matrix estimate: 4 * n^2 / pi^2
    this->cond_est = 4.0 * n * n / std::pow( M_PI, 2 );
  }
};

TYPED_TEST_SUITE_P( SolverTestLaplace2D );

TYPED_TEST_P( SolverTestLaplace2D, DirectSerial )
{
  LinearSolverParameters params = params_DirectSerial();
  params.isSymmetric = true;
  this->test( params );
}

#if !defined(GEOS_USE_CUDA) && !defined(GEOS_USE_HIP)
TYPED_TEST_P( SolverTestLaplace2D, DirectParallel )
{
  this->test( params_DirectParallel() );
}
#endif

TYPED_TEST_P( SolverTestLaplace2D, GMRES_ILU )
{
  this->test( params_GMRES_ILU() );
}

TYPED_TEST_P( SolverTestLaplace2D, CG_SGS )
{
  this->test( params_CG_SGS() );
}

TYPED_TEST_P( SolverTestLaplace2D, CG_AMG )
{
  this->test( params_CG_AMG() );
}

#if defined(GEOS_USE_CUDA) || defined(GEOS_USE_HIP)
REGISTER_TYPED_TEST_SUITE_P( SolverTestLaplace2D,
                             DirectSerial,
                             GMRES_ILU,
                             CG_SGS,
                             CG_AMG );
#else
REGISTER_TYPED_TEST_SUITE_P( SolverTestLaplace2D,
                             DirectSerial,
                             DirectParallel,
                             GMRES_ILU,
                             CG_SGS,
                             CG_AMG );
#endif

#ifdef GEOS_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, SolverTestLaplace2D, TrilinosInterface, );
#endif

#ifdef GEOS_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, SolverTestLaplace2D, HypreInterface, );
#endif

#ifdef GEOS_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, SolverTestLaplace2D, PetscInterface, );
#endif

///////////////////////////////////////////////////////////////////////////////////////

template< typename LAI >
class SolverTestElasticity2D : public SolverTestBase< LAI >
{
public:

  using Base = SolverTestBase< LAI >;
  using Matrix = typename Base::Matrix;
  using Vector = typename Base::Vector;

protected:

  void SetUp() override
  {
    globalIndex constexpr n = 100;
    geos::testing::compute2DElasticityOperator( MPI_COMM_GEOSX, 1.0, 1.0, n, n, 10000., 0.2, this->matrix );
    this->cond_est = 1e4; // not a true condition number estimate, but enough to pass tests
  }
};

TYPED_TEST_SUITE_P( SolverTestElasticity2D );

TYPED_TEST_P( SolverTestElasticity2D, DirectSerial )
{
  this->test( params_DirectSerial() );
}

#if !defined(GEOS_USE_CUDA) && !defined(GEOS_USE_HIP)
TYPED_TEST_P( SolverTestElasticity2D, DirectParallel )
{
  this->test( params_DirectParallel() );
}
#endif

TYPED_TEST_P( SolverTestElasticity2D, GMRES_AMG )
{
  LinearSolverParameters params = params_GMRES_AMG();
  params.amg.separateComponents = true;
  params.dofsPerNode = 2;
  this->test( params );
}

#if defined(GEOS_USE_CUDA) || defined(GEOS_USE_HIP)
REGISTER_TYPED_TEST_SUITE_P( SolverTestElasticity2D,
                             DirectSerial,
                             GMRES_AMG );
#else
REGISTER_TYPED_TEST_SUITE_P( SolverTestElasticity2D,
                             DirectSerial,
                             DirectParallel,
                             GMRES_AMG );
#endif

#ifdef GEOS_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, SolverTestElasticity2D, TrilinosInterface, );
#endif

#ifdef GEOS_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, SolverTestElasticity2D, HypreInterface, );
#endif

#ifdef GEOS_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, SolverTestElasticity2D, PetscInterface, );
#endif

int main( int argc, char * * argv )
{
  geos::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
