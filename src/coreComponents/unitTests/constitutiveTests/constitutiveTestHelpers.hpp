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

#pragma once

// Source includes
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilitySelector.hpp"
#include "constitutive/capillaryPressure/capillaryPressureSelector.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"
#include "unitTests/fluidFlowTests/testFlowUtils.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

namespace geos
{
namespace testing
{

void initializeTable( string const & tableName,
                      array1d< array1d< real64 > > const & coordinates,
                      array1d< real64 > const & values )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  TableFunction & table =
    dynamicCast< TableFunction & >( *functionManager.createChild( TableFunction::catalogName(), tableName ) );
  table.setTableCoordinates( coordinates );
  table.setTableValues( values );
  table.reInitializeFunction();

  table.setInterpolationMethod( TableFunction::InterpolationType::Linear );
}

template< typename MODEL, typename VAR, typename D_VAR_D_SAT >
void testNumericalDerivatives( dataRepository::Group & parent,
                               MODEL & model,
                               arraySlice1d< real64 const > const saturationInput,
                               real64 const perturbParameter,
                               real64 const relTol,
                               string const & varName,
                               VAR && varAccessor,
                               D_VAR_D_SAT && dVar_dSat_accessor )
{
  localIndex const NP = model.numFluidPhases();
  auto const & phases = model.phaseNames();

  // Copy input values into an array with expected layout
  array2d< real64, compflow::LAYOUT_PHASE > saturationValues( 1, NP );
  for( integer i = 0; i < NP; ++i )
  {
    saturationValues[0][i] = saturationInput[i];
  }
  arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const saturation = saturationValues[0];

  // create a clone of the rel perm to run updates on
  std::unique_ptr< constitutive::ConstitutiveBase > modelCopyPtr = model.deliverClone( "fluidCopy", &parent );
  MODEL & modelCopy = dynamicCast< MODEL & >( *modelCopyPtr );

  model.allocateConstitutiveData( model.getParent(), 1 );
  modelCopy.allocateConstitutiveData( model.getParent(), 1 );

  // auto to avoid having to spell out unknown layout permutations
  auto const var = varAccessor( model );
  auto const dPhaseRelPerm_dSat = dVar_dSat_accessor( model );
  auto const varCopy = varAccessor( modelCopy );

  // set the fluid state to current
  constitutive::constitutiveUpdatePassThru( model, [&] ( auto & castedModel )
  {
    typename TYPEOFREF( castedModel ) ::KernelWrapper relPermWrapper = castedModel.createKernelWrapper();
    relPermWrapper.update( 0, 0, saturation );
  } );

  // update saturation and check derivatives
  auto dPhaseRelPerm_dS = testing::invertLayout( dPhaseRelPerm_dSat, NP, NP );

  array2d< real64, compflow::LAYOUT_PHASE > satNew( 1, NP );
  for( integer jp = 0; jp < NP; ++jp )
  {
    real64 const dS = perturbParameter * (saturation[jp] + perturbParameter);
    for( integer ip = 0; ip < NP; ++ip )
    {
      satNew[0][ip] = saturation[ip];
    }
    satNew[0][jp] += dS;

    constitutive::constitutiveUpdatePassThru( modelCopy, [&] ( auto & castedRelPerm )
    {
      typename TYPEOFREF( castedRelPerm ) ::KernelWrapper relPermWrapper = castedRelPerm.createKernelWrapper();
      relPermWrapper.update( 0, 0, satNew[0] );
    } );

    checkDerivative( varCopy.toSliceConst(),
                     var.toSliceConst(),
                     dPhaseRelPerm_dS[jp].toSliceConst(),
                     dS,
                     relTol,
                     varName,
                     "phaseVolFrac[" + phases[jp] + "]",
                     phases );
  }
}

template< typename BASE >
class ConstitutiveTestBase : public ::testing::Test
{
public:
  ConstitutiveTestBase():
    m_node(),
    m_parent( "parent", m_node )
  {
    m_parent.resize( 1 );
  }

  void initialize( BASE & model )
  {
    m_model = &model;
    m_parent.initialize();
    m_parent.initializePostInitialConditions();
  }

protected:
  conduit::Node m_node;
  dataRepository::Group m_parent;
  BASE * m_model;
};

} // namespace testing
} // namespace geos
