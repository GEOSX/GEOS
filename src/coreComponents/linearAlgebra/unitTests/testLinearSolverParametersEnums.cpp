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

#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"

#include <gtest/gtest.h>

using namespace geos;

// The `ENUM_STRING` implementation relies on consistency between the order of the `enum`,
// and the order of the `string` array provided. Since this consistency is not enforced, it can be corrupted anytime.
// This unit test aims at preventing from this implicit relationship to bring a bug.

TEST( LinearSolverParametersEnums, SolverType )
{
  using EnumType = LinearSolverParameters::SolverType;

  ASSERT_EQ( "direct", toString( EnumType::direct ) );
  ASSERT_EQ( "cg", toString( EnumType::cg ) );
  ASSERT_EQ( "gmres", toString( EnumType::gmres ) );
  ASSERT_EQ( "fgmres", toString( EnumType::fgmres ) );
  ASSERT_EQ( "bicgstab", toString( EnumType::bicgstab ) );
  ASSERT_EQ( "preconditioner", toString( EnumType::preconditioner ) );
}


TEST( LinearSolverParametersEnums, PreconditionerType )
{
  using EnumType = LinearSolverParameters::PreconditionerType;

  ASSERT_EQ( "none", toString( EnumType::none ) );
  ASSERT_EQ( "jacobi", toString( EnumType::jacobi ) );
  ASSERT_EQ( "l1jacobi", toString( EnumType::l1jacobi ) );
  ASSERT_EQ( "fgs", toString( EnumType::fgs ) );
  ASSERT_EQ( "sgs", toString( EnumType::sgs ) );
  ASSERT_EQ( "l1sgs", toString( EnumType::l1sgs ) );
  ASSERT_EQ( "chebyshev", toString( EnumType::chebyshev ) );
  ASSERT_EQ( "iluk", toString( EnumType::iluk ) );
  ASSERT_EQ( "ilut", toString( EnumType::ilut ) );
  ASSERT_EQ( "icc", toString( EnumType::ic ) ); // Notice the discrepancy here
  ASSERT_EQ( "ict", toString( EnumType::ict ) );
  ASSERT_EQ( "amg", toString( EnumType::amg ) );
  ASSERT_EQ( "mgr", toString( EnumType::mgr ) );
  ASSERT_EQ( "block", toString( EnumType::block ) );
  ASSERT_EQ( "direct", toString( EnumType::direct ) );
  ASSERT_EQ( "bgs", toString( EnumType::bgs ) );
}


TEST( LinearSolverParametersEnums, DirectColPerm )
{
  using EnumType = LinearSolverParameters::Direct::ColPerm;

  ASSERT_EQ( "none", toString( EnumType::none ) );
  ASSERT_EQ( "MMD_AtplusA", toString( EnumType::MMD_AtplusA ) );
  ASSERT_EQ( "MMD_AtA", toString( EnumType::MMD_AtA ) );
  ASSERT_EQ( "colAMD", toString( EnumType::colAMD ) );
  ASSERT_EQ( "metis", toString( EnumType::metis ) );
  ASSERT_EQ( "parmetis", toString( EnumType::parmetis ) );
}


TEST( LinearSolverParametersEnums, DirectRowPerm )
{
  using EnumType = LinearSolverParameters::Direct::RowPerm;

  ASSERT_EQ( "none", toString( EnumType::none ) );
  ASSERT_EQ( "mc64", toString( EnumType::mc64 ) );
}


TEST( LinearSolverParametersEnums, MGRStrategyType )
{
  using EnumType = LinearSolverParameters::MGR::StrategyType;

  ASSERT_EQ( "invalid", toString( EnumType::invalid ) );
  ASSERT_EQ( "singlePhaseReservoirFVM", toString( EnumType::singlePhaseReservoirFVM ) );
  ASSERT_EQ( "singlePhaseHybridFVM", toString( EnumType::singlePhaseHybridFVM ) );
  ASSERT_EQ( "singlePhaseReservoirHybridFVM", toString( EnumType::singlePhaseReservoirHybridFVM ) );
  ASSERT_EQ( "singlePhasePoromechanics", toString( EnumType::singlePhasePoromechanics ) );
  ASSERT_EQ( "hybridSinglePhasePoromechanics", toString( EnumType::hybridSinglePhasePoromechanics ) );
  ASSERT_EQ( "compositionalMultiphaseFVM", toString( EnumType::compositionalMultiphaseFVM ) );
  ASSERT_EQ( "compositionalMultiphaseHybridFVM", toString( EnumType::compositionalMultiphaseHybridFVM ) );
  ASSERT_EQ( "compositionalMultiphaseReservoirFVM", toString( EnumType::compositionalMultiphaseReservoirFVM ) );
  ASSERT_EQ( "compositionalMultiphaseReservoirHybridFVM", toString( EnumType::compositionalMultiphaseReservoirHybridFVM ) );
  ASSERT_EQ( "multiphasePoromechanics", toString( EnumType::multiphasePoromechanics ) );
  ASSERT_EQ( "hydrofracture", toString( EnumType::hydrofracture ) );
  ASSERT_EQ( "lagrangianContactMechanics", toString( EnumType::lagrangianContactMechanics ) );
}


TEST( LinearSolverParametersEnums, AMGCycleType )
{
  using EnumType = LinearSolverParameters::AMG::CycleType;

  ASSERT_EQ( "V", toString( EnumType::V ) );
  ASSERT_EQ( "W", toString( EnumType::W ) );
}


TEST( LinearSolverParametersEnums, AMGPreOrPost )
{
  using EnumType = LinearSolverParameters::AMG::PreOrPost;

  ASSERT_EQ( "pre", toString( EnumType::pre ) );
  ASSERT_EQ( "post", toString( EnumType::post ) );
  ASSERT_EQ( "both", toString( EnumType::both ) );
}


TEST( LinearSolverParametersEnums, AMGSmootherType )
{
  using EnumType = LinearSolverParameters::AMG::SmootherType;

  ASSERT_EQ( "default", toString( EnumType::default_ ) );
  ASSERT_EQ( "jacobi", toString( EnumType::jacobi ) );
  ASSERT_EQ( "l1jacobi", toString( EnumType::l1jacobi ) );
  ASSERT_EQ( "fgs", toString( EnumType::fgs ) );
  ASSERT_EQ( "bgs", toString( EnumType::bgs ) );
  ASSERT_EQ( "sgs", toString( EnumType::sgs ) );
  ASSERT_EQ( "l1sgs", toString( EnumType::l1sgs ) );
  ASSERT_EQ( "chebyshev", toString( EnumType::chebyshev ) );
  ASSERT_EQ( "ilu0", toString( EnumType::ilu0 ) );
  ASSERT_EQ( "ilut", toString( EnumType::ilut ) );
  ASSERT_EQ( "ic0", toString( EnumType::ic0 ) );
  ASSERT_EQ( "ict", toString( EnumType::ict ) );
}


TEST( LinearSolverParametersEnums, AMGNullSpaceType )
{
  using EnumType = LinearSolverParameters::AMG::NullSpaceType;

  ASSERT_EQ( "constantModes", toString( EnumType::constantModes ) );
  ASSERT_EQ( "rigidBodyModes", toString( EnumType::rigidBodyModes ) );
}

int main( int argc, char * * argv )
{
  geos::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
