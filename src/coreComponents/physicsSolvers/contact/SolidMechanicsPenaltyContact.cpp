/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsPenaltyContact.cpp
 *
 */

#include "SolidMechanicsPenaltyContact.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/FrictionSelector.hpp"


#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;
using namespace finiteElement;

SolidMechanicsPenaltyContact::SolidMechanicsPenaltyContact( const string & name,
                                                            Group * const parent ):
  ContactSolverBase( name, parent )
{
  registerWrapper( viewKeyStruct::contactPenaltyStiffnessString(), &m_contactPenaltyStiffness ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Value of the penetration penalty stiffness. Units of Pressure/length" );
}

SolidMechanicsPenaltyContact::~SolidMechanicsPenaltyContact()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanicsPenaltyContact::setupSystem( DomainPartition & domain,
                                                DofManager & dofManager,
                                                CRSMatrix< real64, globalIndex > & localMatrix,
                                                ParallelVector & rhs,
                                                ParallelVector & solution,
                                                bool const setSparsity )
{
  GEOS_MARK_FUNCTION;
  PhysicsSolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, false );

  SparsityPattern< globalIndex > sparsityPattern( dofManager.numLocalDofs(),
                                                  dofManager.numGlobalDofs(),
                                                  8*8*3*1.2 );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    arrayView1d< globalIndex const > const
    dofNumber = nodeManager.getReference< globalIndex_array >( dofManager.getKey( solidMechanics::totalDisplacement::key() ) );


    ElementRegionManager const & elemManager = mesh.getElemManager();
    array1d< string > allFaceElementRegions;
    elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & elemRegion )
    {
      allFaceElementRegions.emplace_back( elemRegion.getName() );
    } );

    finiteElement::
      fillSparsity< FaceElementSubRegion,
                    solidMechanicsLagrangianFEMKernels::ImplicitSmallStrainQuasiStatic >( mesh,
                                                                                          allFaceElementRegions,
                                                                                          this->getDiscretizationName(),
                                                                                          dofNumber,
                                                                                          dofManager.rankOffset(),
                                                                                          sparsityPattern );

    finiteElement::fillSparsity< CellElementSubRegion,
                                 solidMechanicsLagrangianFEMKernels::ImplicitSmallStrainQuasiStatic >( mesh,
                                                                                                       regionNames,
                                                                                                       this->getDiscretizationName(),
                                                                                                       dofNumber,
                                                                                                       dofManager.rankOffset(),
                                                                                                       sparsityPattern );


  } );

  sparsityPattern.compress();
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( sparsityPattern ) );
}

void SolidMechanicsPenaltyContact::assembleSystem( real64 const time,
                                                   real64 const dt,
                                                   DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  synchronizeFractureState( domain );

  SolidMechanicsLagrangianFEM::assembleSystem( time,
                                               dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );

  assembleContact( domain, dofManager, localMatrix, localRhs );
}

void SolidMechanicsPenaltyContact::assembleContact( DomainPartition & domain,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager & nodeManager = mesh.getNodeManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    solidMechanics::arrayViewConst2dLayoutTotalDisplacement const u =
      nodeManager.getField< solidMechanics::totalDisplacement >();
    arrayView2d< real64 > const fc = nodeManager.getField< solidMechanics::contactForce >();
    fc.zero();

    arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    string const dofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    arrayView1d< globalIndex > const nodeDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );
    globalIndex const rankOffset = dofManager.rankOffset();

    // TODO: this bound may need to change
    constexpr localIndex maxNodexPerFace = 4;
    constexpr localIndex maxDofPerElem = maxNodexPerFace * 3 * 2;

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      real64 const contactStiffness = m_contactPenaltyStiffness;

      arrayView1d< real64 > const area = subRegion.getElementArea();
      arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

      // TODO: use parallel policy?
      forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        localIndex const kf0 = elemsToFaces[kfe][0], kf1 = elemsToFaces[kfe][1];
        real64 Nbar[ 3 ] = { faceNormal[kf0][0] - faceNormal[kf1][0],
                             faceNormal[kf0][1] - faceNormal[kf1][1],
                             faceNormal[kf0][2] - faceNormal[kf1][2] };

        LvArray::tensorOps::normalize< 3 >( Nbar );

        localIndex const numNodesPerFace=facesToNodes.sizeOfArray( kf0 );
        real64 const Ja = area[kfe] / numNodesPerFace;

        stackArray1d< globalIndex, maxDofPerElem > rowDOF( numNodesPerFace*3*2 );
        stackArray1d< real64, maxDofPerElem > nodeRHS( numNodesPerFace*3*2 );
        stackArray2d< real64, maxDofPerElem *maxDofPerElem > dRdP( numNodesPerFace*3*2, numNodesPerFace*3*2 );

        for( localIndex a=0; a<numNodesPerFace; ++a )
        {
          real64 penaltyForce[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( Nbar );
          localIndex const node0 = facesToNodes[kf0][a];
          localIndex const node1 = facesToNodes[kf1][ a==0 ? a : numNodesPerFace-a ];
          real64 gap[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( u[node1] );
          LvArray::tensorOps::subtract< 3 >( gap, u[node0] );
          real64 const gapNormal = LvArray::tensorOps::AiBi< 3 >( gap, Nbar );

          for( int i=0; i<3; ++i )
          {
            rowDOF[3*a+i]                     = nodeDofNumber[node0]+i;
            rowDOF[3*(numNodesPerFace + a)+i] = nodeDofNumber[node1]+i;
          }

          if( gapNormal < 0 )
          {
            LvArray::tensorOps::scale< 3 >( penaltyForce, -contactStiffness * gapNormal * Ja );
            for( int i=0; i<3; ++i )
            {
              LvArray::tensorOps::subtract< 3 >( fc[node0], penaltyForce );
              LvArray::tensorOps::add< 3 >( fc[node1], penaltyForce );
              nodeRHS[3*a+i]                     -= penaltyForce[i];
              nodeRHS[3*(numNodesPerFace + a)+i] += penaltyForce[i];

              dRdP( 3*a+i, 3*a+i )                                         -= contactStiffness * Ja * Nbar[i] * Nbar[i];
              dRdP( 3*a+i, 3*(numNodesPerFace + a)+i )                     += contactStiffness * Ja * Nbar[i] * Nbar[i];
              dRdP( 3*(numNodesPerFace + a)+i, 3*a+i )                     += contactStiffness * Ja * Nbar[i] * Nbar[i];
              dRdP( 3*(numNodesPerFace + a)+i, 3*(numNodesPerFace + a)+i ) -= contactStiffness * Ja * Nbar[i] * Nbar[i];
            }
          }
        }

        for( localIndex idof = 0; idof < numNodesPerFace*3*2; ++idof )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( rowDOF[idof] - rankOffset );

          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                      rowDOF.data(),
                                                                      dRdP[idof].dataIfContiguous(),
                                                                      numNodesPerFace*3*2 );
            RAJA::atomicAdd( serialAtomic{}, &localRhs[localRow], nodeRHS[idof] );
          }
        }
      } );
    } );
  } );
}


REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SolidMechanicsPenaltyContact, string const &, Group * const )

} /* namespace geos */
