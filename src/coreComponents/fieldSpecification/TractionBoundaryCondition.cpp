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
 * @file TractionBoundaryCondition.cpp
 */

#include "TractionBoundaryCondition.hpp"

#include "functions/TableFunction.hpp"

namespace geos
{

using namespace dataRepository;

TractionBoundaryCondition::TractionBoundaryCondition( string const & name, Group * parent ):
  FieldSpecificationBase( name, parent ),
  m_tractionType( TractionType::vector ),
  m_inputStress{}//,
//  m_stressFunctionNames(),
//  m_useStressFunctions(false),
//  m_stressFunctions{nullptr}
{
  registerWrapper( viewKeyStruct::tractionTypeString(), &m_tractionType ).
    setApplyDefaultValue( TractionType::vector ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Type of traction boundary condition. Options are:\n" +
                    toString( TractionType::vector ) + " - traction is applied to the faces as specified from the scale and direction,\n" +
                    toString( TractionType::normal ) + " - traction is applied to the faces as a pressure specified from the product of scale and the outward face normal,\n" +
                    toString( TractionType::stress ) + " - traction is applied to the faces as specified by the inner product of input stress and face normal." );

  registerWrapper( viewKeyStruct::inputStressString(), &m_inputStress ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( string( "Input stress for " ) + viewKeyStruct::tractionTypeString() + " = " + toString( TractionType::stress ) );

//  registerWrapper( viewKeyStruct::stressFunctionString(), &m_stressFunctionNames ).
//    setInputFlag( InputFlags::OPTIONAL ).
//    setDescription( string("Function names for description of stress for ") + viewKeyStruct::tractionTypeString() +
//                    " = " + toString( TractionType::stress ) + ". Overrides " + viewKeyStruct::inputStressString() + "." );

  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).
    setInputFlag( InputFlags::FALSE );
  setFieldName( catalogName() );

  getWrapper< int >( FieldSpecificationBase::viewKeyStruct::componentString() ).
    setInputFlag( InputFlags::FALSE );
}


void TractionBoundaryCondition::postInputInitialization()
{
  if( m_tractionType == TractionType::vector )
  {
    GEOS_ERROR_IF( LvArray::tensorOps::l2Norm< 3 >( getDirection() ) < 1e-20,
                   getDataContext() << ": " << viewKeyStruct::directionString() << " is required for " <<
                   viewKeyStruct::tractionTypeString() << " = " << TractionType::vector <<
                   ", but appears to be unspecified" );
  }
  else
  {
    GEOS_LOG_RANK_0_IF( LvArray::tensorOps::l2Norm< 3 >( getDirection() ) > 1e-20,
                        getDataContext() << ": " << viewKeyStruct::directionString() << " is not required unless " <<
                        viewKeyStruct::tractionTypeString() << " = " << TractionType::vector <<
                        ", but appears to be specified" );
  }

  bool const inputStressRead = getWrapper< R2SymTensor >( viewKeyStruct::inputStressString() ).getSuccessfulReadFromInput();

  GEOS_LOG_RANK_0_IF( inputStressRead && m_tractionType != TractionType::stress,
                      getDataContext() << ": " << viewKeyStruct::inputStressString() << " is specified, but " <<
                      viewKeyStruct::tractionTypeString() << " != " << TractionType::stress <<
                      ", so value of " << viewKeyStruct::inputStressString() << " is unused." );

  GEOS_ERROR_IF( !inputStressRead && m_tractionType == TractionType::stress,
                 getDataContext() << ": " << viewKeyStruct::tractionTypeString() << " = " << TractionType::stress <<
                 ", but " << viewKeyStruct::inputStressString() << " is not specified." );


//  localIndex const numStressFunctionsNames = m_stressFunctionNames.size();
//  GEOS_ERROR_IF( numStressFunctionsNames > 0 && numStressFunctionsNames<6,
//                  "Either 0 or 6 stress functions must be specified using stressFunctions" );
//
//  if( numStressFunctionsNames==6 )
//  {
//    m_useStressFunctions = true;
//  }
}

void TractionBoundaryCondition::initializePreSubGroups()
{
//  FunctionManager const & functionManager = getGlobalState().getFunctionManager();
//
//  if( m_useStressFunctions )
//  {
//    for( localIndex i=0; i<6; ++i )
//    {
//      m_stressFunctions[i] = &( functionManager.getGroup< TableFunction >( m_stressFunctionNames[i] ) );
//    }
//  }
}


void TractionBoundaryCondition::launch( real64 const time,
                                        arrayView1d< globalIndex const > const blockLocalDofNumber,
                                        globalIndex const dofRankOffset,
                                        FaceManager const & faceManager,
                                        SortedArrayView< localIndex const > const & targetSet,
                                        arrayView1d< real64 > const & localRhs ) const
{


  arrayView1d< real64 const > const faceArea  = faceManager.faceArea();
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  FunctionManager const & functionManager = FunctionManager::getInstance();

  string const & functionName = this->getFunctionName();

  globalIndex_array nodeDOF;
  real64_array nodeRHS;

  R1Tensor const direction = this->getDirection();
  TractionType const tractionType = m_tractionType;
  R2SymTensor inputStress = m_inputStress;

  real64 tractionMagnitude0;
  array1d< real64 > tractionMagnitudeArray( targetSet.size() );
  bool spatialFunction = false;

  if( functionName.empty() )
  {
    tractionMagnitude0 = this->getScale();
  }
  else
  {
    FunctionBase const & function = functionManager.getGroup< FunctionBase >( functionName );
    if( function.isFunctionOfTime() == 2 )
    {
      tractionMagnitude0 = this->getScale() * function.evaluate( &time );
    }
    else
    {
      tractionMagnitudeArray.setName( getName() + " function results" );
      function.evaluate( faceManager, time, targetSet, tractionMagnitudeArray );
      spatialFunction = true;
      tractionMagnitude0 = 1e99;
    }
  }

  arrayView1d< real64 const > const tractionMagnitudeArrayView = tractionMagnitudeArray;

  {
    forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localIndex const kf = targetSet[ i ];
      localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );

      // TODO consider dispatch if appropriate
      real64 const tractionMagnitude = spatialFunction ? tractionMagnitudeArrayView[i] : tractionMagnitude0;

      real64 traction[3] = { 0 };
      // TODO consider dispatch if appropriate
      if( tractionType == TractionType::vector )
      {
        traction[0] = tractionMagnitude * direction[0];
        traction[1] = tractionMagnitude * direction[1];
        traction[2] = tractionMagnitude * direction[2];
      }
      else
      {
        real64 const temp[3] = { tractionMagnitude * faceNormal( kf, 0 ),
                                 tractionMagnitude * faceNormal( kf, 1 ),
                                 tractionMagnitude * faceNormal( kf, 2 ) };

        if( tractionType == TractionType::normal )
        {
          traction[0] = temp[0];
          traction[1] = temp[1];
          traction[2] = temp[2];
        }
        else if( tractionType == TractionType::stress )
        {
          traction[0] = inputStress[0] * temp[0] + inputStress[5] * temp[1] + inputStress[4] * temp[2];
          traction[1] = inputStress[5] * temp[0] + inputStress[1] * temp[1] + inputStress[3] * temp[2];
          traction[2] = inputStress[4] * temp[0] + inputStress[3] * temp[1] + inputStress[2] * temp[2];
        }

      }

      // TODO replace with proper FEM integration
      traction[0] *= faceArea[kf] / numNodes;
      traction[1] *= faceArea[kf] / numNodes;
      traction[2] *= faceArea[kf] / numNodes;

      for( localIndex a=0; a<numNodes; ++a )
      {
        localIndex const dof = blockLocalDofNumber[ faceToNodeMap( kf, a ) ] - dofRankOffset;
        if( dof < 0 || dof >= localRhs.size() )
          continue;
        RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof+0], traction[0] );
        RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof+1], traction[1] );
        RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof+2], traction[2] );
      }
    } );
  }
}

REGISTER_CATALOG_ENTRY( FieldSpecificationBase, TractionBoundaryCondition, string const &, Group * const )


} /* namespace geos */
