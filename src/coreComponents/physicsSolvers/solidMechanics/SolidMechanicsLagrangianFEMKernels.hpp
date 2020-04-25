/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsLagrangianFEMKernels.hpp
 */

#pragma once

#include "../PhysicsLoopInterface.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteElement/ElementLibrary/FiniteElementBase.h"
#include "finiteElement/Kinematics.h"
#include "TimeIntegrationOption.hpp"

namespace geosx
{

namespace SolidMechanicsLagrangianFEMKernels
{

inline void velocityUpdate( arrayView2d< real64, nodes::ACCELERATION_USD > const & acceleration,
                            arrayView2d< real64, nodes::VELOCITY_USD > const & velocity,
                            real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  localIndex const N = acceleration.size( 0 );
  forAll< parallelDevicePolicy<> >( N, [=] GEOSX_DEVICE ( localIndex const i )
  {
    for( int j = 0; j < 3; ++j )
    {
      velocity( i, j ) += dt * acceleration( i, j );
      acceleration( i, j ) = 0;
    }
  } );
}

inline void velocityUpdate( arrayView2d< real64, nodes::ACCELERATION_USD > const & acceleration,
                            arrayView1d< real64 const > const & mass,
                            arrayView2d< real64, nodes::VELOCITY_USD > const & velocity,
                            real64 const dt,
                            SortedArrayView< localIndex const > const & indices )
{
  GEOSX_MARK_FUNCTION;

  forAll< parallelDevicePolicy<> >( indices.size(), [=] GEOSX_DEVICE ( localIndex const i )
  {
    localIndex const a = indices[ i ];
    for( int j = 0; j < 3; ++j )
    {
      acceleration( a, j ) /= mass[ a ];
      velocity( a, j ) += dt * acceleration( a, j );
    }
  } );
}

inline void displacementUpdate( arrayView2d< real64 const, nodes::VELOCITY_USD > const & velocity,
                                arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & uhat,
                                arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & u,
                                real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  localIndex const N = velocity.size( 0 );
  forAll< parallelDevicePolicy<> >( N, [=] GEOSX_DEVICE ( localIndex const i )
  {
    for( int j = 0; j < 3; ++j )
    {
      uhat( i, j ) = velocity( i, j ) * dt;
      u( i, j ) += uhat( i, j );
    }
  } );
}

template< int N >
inline void Integrate( const R2SymTensor & fieldvar,
                       arraySlice1d< R1Tensor const > const & dNdX,
                       real64 const & detJ,
                       real64 const & detF,
                       const R2Tensor & fInv,
                       R1Tensor * GEOSX_RESTRICT const result )
{
  real64 const integrationFactor = detJ * detF;

  R2Tensor P;
  P.AijBkj( fieldvar, fInv );
  P *= integrationFactor;

  for( int a=0; a<N; ++a )    // loop through all shape functions in element
  {
    result[a].minusAijBj( P, dNdX[a] );
  }
}

/**
 * @brief Function to select which templated kernel function to call.
 * @tparam KERNELWRAPPER A struct or class that contains the following method
 *  "Launch<NUM_NODES_PER_ELEM, NUM_QUADRATURE_POINTS, CONSTITUTIVE_TYPE>( CONSTITUTIVE_TYPE *, PARAMS... )"
 * @tparam PARAMS Variadic parameter pack to pass arguments to Launch function.
 * @param NUM_NODES_PER_ELEM The number of nodes in an element.
 * @param NUM_QUADRATURE_POINTS The number of quadrature points in an element.
 * @param params Variadic parameter list to hold all parameters that are forwarded to the kernel function.
 * @return Depends on the kernel.
 */
template< typename KERNELWRAPPER, typename ... PARAMS >
inline real64
ElementKernelLaunchSelector( localIndex NUM_NODES_PER_ELEM,
                             localIndex NUM_QUADRATURE_POINTS,
                             constitutive::ConstitutiveBase * const constitutiveRelation,
                             PARAMS && ... params )
{
  real64 rval = 0;

  using namespace constitutive;

  ConstitutivePassThru< SolidBase >::Execute( constitutiveRelation,
                                              [&]( auto * const constitutive )
  {
    using CONSTITUTIVE_TYPE = TYPEOFPTR( constitutive );
    if( NUM_NODES_PER_ELEM==8 && NUM_QUADRATURE_POINTS==8 )
    {
      rval = KERNELWRAPPER::template Launch< 8, 8, CONSTITUTIVE_TYPE >( constitutive,
                                                                        std::forward< PARAMS >( params )... );
    }
    else if( NUM_NODES_PER_ELEM==4 && NUM_QUADRATURE_POINTS==1 )
    {
      rval = KERNELWRAPPER::template Launch< 4, 1, CONSTITUTIVE_TYPE >( constitutive,
                                                                        std::forward< PARAMS >( params )... );
    }
  } );
  return rval;
}

/**
 * @struct Structure to wrap templated function that implements the explicit time integration kernel.
 */
struct ExplicitKernel
{
  /**
   * @brief Launch of the element processing kernel for explicit time integration.
   * @tparam NUM_NODES_PER_ELEM The number of nodes/dof per element.
   * @tparam NUM_QUADRATURE_POINTS The number of quadrature points per element.
   * @tparam CONSTITUTIVE_TYPE the type of the constitutive relation that is being used.
   * @param A pointer to the constitutive relation that is being used.
   * @param elementList The list of elements to be processed
   * @param elemsToNodes The map from the elements to the nodes that form that element.
   * @param dNdX The derivatives of the shape functions wrt the reference configuration.
   * @param detJ The determinant of the transformation matrix (Jacobian) to the parent element.
   * @param u The nodal array of total displacements.
   * @param vel The nodal array of velocity.
   * @param acc The nodal array of force/acceleration.
   * @param stress The stress at each element quadrature point
   * @param dt The timestep
   * @return The achieved timestep.
   */
  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS, typename CONSTITUTIVE_TYPE >
  static inline real64
  Launch( CONSTITUTIVE_TYPE * const constitutiveRelation,
          SortedArrayView< localIndex const > const & elementList,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView3d< R1Tensor const > const & dNdX,
          arrayView2d< real64 const > const & detJ,
          arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u,
          arrayView2d< real64 const, nodes::VELOCITY_USD > const & vel,
          arrayView2d< real64, nodes::ACCELERATION_USD > const & acc,
          real64 const dt )
  {

    typename CONSTITUTIVE_TYPE::KernelWrapper const & constitutive = constitutiveRelation->createKernelWrapper();

    forAll< serialPolicy >( elementList.size(), [=] ( localIndex const i )
    {
      localIndex const k = elementList[ i ];
      R1Tensor v_local[NUM_NODES_PER_ELEM];
      R1Tensor u_local[NUM_NODES_PER_ELEM];
      R1Tensor f_local[NUM_NODES_PER_ELEM];

      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        u_local[ a ] = u[ nodeIndex ];
        v_local[ a ] = vel[ nodeIndex ];
      }

      //Compute Quadrature
      for( localIndex q = 0; q<NUM_QUADRATURE_POINTS; ++q )
      {
        R2Tensor dUhatdX, dUdX;
        CalculateGradients< NUM_NODES_PER_ELEM >( dUhatdX, dUdX, v_local, u_local, dNdX[k][q] );
        dUhatdX *= dt;

        R2Tensor F, Ldt, fInv;

        // calculate du/dX
        F = dUhatdX;
        F *= 0.5;
        F += dUdX;
        F.PlusIdentity( 1.0 );
        fInv.Inverse( F );

        // chain rule: calculate dv/du = dv/dX * dX/du
        Ldt.AijBjk( dUhatdX, fInv );

        // calculate gradient (end of step)
        F = dUhatdX;
        F += dUdX;
        F.PlusIdentity( 1.0 );
        real64 detF = F.Det();
        fInv.Inverse( F );


        R2Tensor Rot;
        R2SymTensor Dadt;
        HughesWinget( Rot, Dadt, Ldt );

        constitutive.HypoElastic( k, q, Dadt.Data(), Rot );

        Integrate< NUM_NODES_PER_ELEM >( constitutive.m_stress[k][q], dNdX[k][q], detJ[k][q], detF, fInv, f_local );
      }//quadrature loop

      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        RAJA::atomicAdd< serialAtomic >( &acc( nodeIndex, 0 ), f_local[ a ][ 0 ] );
        RAJA::atomicAdd< serialAtomic >( &acc( nodeIndex, 1 ), f_local[ a ][ 1 ] );
        RAJA::atomicAdd< serialAtomic >( &acc( nodeIndex, 2 ), f_local[ a ][ 2 ] );
      }

    } );

    return dt;
  }


  static inline real64
  CalculateSingleNodalForce( localIndex const k,
                             localIndex const targetNode,
                             localIndex const numQuadraturePoints,
                             arrayView3d< R1Tensor const > const & dNdX,
                             arrayView2d< real64 const > const & detJ,
                             arrayView3d< real64 const, solid::STRESS_USD > const & stress,
                             R1Tensor & force )
  {
    GEOSX_MARK_FUNCTION;
    localIndex const & a = targetNode;

    //Compute Quadrature
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      force[ 0 ] -= ( stress( k, q, 0 ) * dNdX( k, q, a )[ 0 ] +
                      stress( k, q, 5 ) * dNdX( k, q, a )[ 1 ] +
                      stress( k, q, 4 ) * dNdX( k, q, a )[ 2 ] ) * detJ( k, q );
      force[ 1 ] -= ( stress( k, q, 5 ) * dNdX( k, q, a )[ 0 ] +
                      stress( k, q, 1 ) * dNdX( k, q, a )[ 1 ] +
                      stress( k, q, 3 ) * dNdX( k, q, a )[ 2 ] ) * detJ( k, q );
      force[ 2 ] -= ( stress( k, q, 4 ) * dNdX( k, q, a )[ 0 ] +
                      stress( k, q, 3 ) * dNdX( k, q, a )[ 1 ] +
                      stress( k, q, 2 ) * dNdX( k, q, a )[ 2 ] ) * detJ( k, q );

    }//quadrature loop

    return 0;
  }

};

/**
 * @struct Structure to wrap templated function that implements the implicit time integration kernel.
 */
struct ImplicitKernel
{

  /**
   * @brief Launch of the element processing kernel for implicit time integration.
   * @tparam NUM_NODES_PER_ELEM The number of nodes/dof per element.
   * @tparam NUM_QUADRATURE_POINTS The number of quadrature points per element.
   * @tparam CONSTITUTIVE_TYPE the type of the constitutive relation that is being used.
   * @param constitutiveRelation A pointer to the constitutive relation that is being used.
   * @param numElems The number of elements the kernel will process.
   * @param dt The timestep.
   * @param dNdX The derivatives of the shape functions wrt the reference configuration.
   * @param detJ The determinant of the transformation matrix (Jacobian) to the parent element.
   * @param fe A pointer to the finite element class used in this kernel.
   * @param elemGhostRank An array containing the values of the owning ranks for ghost elements.
   * @param elemsToNodes The map from the elements to the nodes that form that element.
   * @param globalDofNumber The map from localIndex to the globalDOF number.
   * @param disp The array of total displacements.
   * @param uhat The array of incremental displacements (displacement for this step).
   * @param vtilde The array for the velocity predictor.
   * @param uhattilde The array for the incremental displacement predictor.
   * @param density The array containing the density
   * @param fluidPressure Array containing element fluid pressure at the beginning of the step.
   * @param deltaFluidPressure Array containing the change in element fluid pressure over this step.
   * @param biotCoefficient The biotCoefficient used to calculate effective stress.
   * @param tiOption The time integration option used for the integration.
   * @param stiffnessDamping The stiffness damping coefficient for the Newmark method assuming Rayleigh damping.
   * @param massDamping The mass damping coefficient for the Newmark method assuming Rayleigh damping.
   * @param newmarkBeta The value of \beta in the Newmark update.
   * @param newmarkGamma The value of \gamma in the Newmark update.
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix sparse matrix containing the derivatives of the residual wrt displacement
   * @param rhs parallel vector containing the global residual
   * @return The maximum nodal force contribution from all elements.
   */
  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS, typename CONSTITUTIVE_TYPE >
  static inline real64
  Launch( CONSTITUTIVE_TYPE * const GEOSX_UNUSED_PARAM( constitutiveRelation ),
          localIndex const GEOSX_UNUSED_PARAM( numElems ),
          real64 const GEOSX_UNUSED_PARAM( dt ),
          arrayView3d< R1Tensor const > const & GEOSX_UNUSED_PARAM( dNdX ),
          arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM( detJ ),
          FiniteElementBase const * const GEOSX_UNUSED_PARAM( fe ),
          arrayView1d< integer const > const & GEOSX_UNUSED_PARAM( elemGhostRank ),
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & GEOSX_UNUSED_PARAM( elemsToNodes ),
          arrayView1d< globalIndex const > const & GEOSX_UNUSED_PARAM( globalDofNumber ),
          arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & GEOSX_UNUSED_PARAM( disp ),
          arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & GEOSX_UNUSED_PARAM( uhat ),
          arrayView1d< R1Tensor const > const & GEOSX_UNUSED_PARAM( vtilde ),
          arrayView1d< R1Tensor const > const & GEOSX_UNUSED_PARAM( uhattilde ),
          arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM( density ),
          arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( fluidPressure ),
          arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( deltaFluidPressure ),
          real64 const GEOSX_UNUSED_PARAM( biotCoefficient ),
          TimeIntegrationOption const GEOSX_UNUSED_PARAM( tiOption ),
          real64 const GEOSX_UNUSED_PARAM( stiffnessDamping ),
          real64 const GEOSX_UNUSED_PARAM( massDamping ),
          real64 const GEOSX_UNUSED_PARAM( newmarkBeta ),
          real64 const GEOSX_UNUSED_PARAM( newmarkGamma ),
          R1Tensor const & GEOSX_UNUSED_PARAM( gravityVector ),
          DofManager const * const GEOSX_UNUSED_PARAM( dofManager ),
          ParallelMatrix * const GEOSX_UNUSED_PARAM( matrix ),
          ParallelVector * const GEOSX_UNUSED_PARAM( rhs ) )
  {
    GEOSX_ERROR( "SolidMechanicsLagrangianFEM::ImplicitElementKernelWrapper::Launch() not implemented" );
    return 0;
  }

};

class QuasiStatic
{
public:
  using Base = physicsLoopInterface::FiniteElementRegionLoop;
  static constexpr int numTestDofPerSP = 3;
  static constexpr int numTrialDofPerSP = 3;

  struct Parameters : public Base::Parameters
  {
    Parameters( real64 const inputGravityVector[3] ):
      m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] }
    {}

    real64 const m_gravityVector[3];
  };


  template< int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
            int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM >
  struct StackVariables : Base::StackVariables< NUM_TEST_SUPPORT_POINTS_PER_ELEM*numTestDofPerSP,
                                                NUM_TRIAL_SUPPORT_POINTS_PER_ELEM*numTrialDofPerSP >
  {
public:
    using StackVariablesBase = Base::StackVariables< NUM_TEST_SUPPORT_POINTS_PER_ELEM*numTestDofPerSP,
                                                     NUM_TRIAL_SUPPORT_POINTS_PER_ELEM*numTrialDofPerSP >;
    using StackVariablesBase::numRows;
    using StackVariablesBase::numCols;
    static constexpr int numNodes = NUM_TEST_SUPPORT_POINTS_PER_ELEM;


//      GEOSX_HOST_DEVICE
    StackVariables():
      StackVariablesBase(),
      u_local{ { 0.0, 0.0, 0.0} },
      uhat_local{ { 0.0, 0.0, 0.0} },
      constitutiveStiffness{ {0.0} }
    {}

    R1Tensor u_local[numNodes];
    R1Tensor uhat_local[numNodes];
    real64 constitutiveStiffness[6][6];
  };

  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_NODES_PER_ELEM,
            int >
  using SparsityKernels = Base::Kernels< SUBREGION_TYPE,
                                         CONSTITUTIVE_TYPE,
                                         NUM_NODES_PER_ELEM,
                                         NUM_NODES_PER_ELEM,
                                         numTestDofPerSP,
                                         numTrialDofPerSP >
  ;

  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_NODES_PER_ELEM,
            int >
  class Kernels : public Base::Kernels< SUBREGION_TYPE,
                                        CONSTITUTIVE_TYPE,
                                        NUM_NODES_PER_ELEM,
                                        NUM_NODES_PER_ELEM,
                                        numTestDofPerSP,
                                        numTrialDofPerSP >
  {
public:
    using KernelBase = Base::Kernels< SUBREGION_TYPE,
                                      CONSTITUTIVE_TYPE,
                                      NUM_NODES_PER_ELEM,
                                      NUM_NODES_PER_ELEM,
                                      numTestDofPerSP,
                                      numTrialDofPerSP >;

    static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;

    using KernelBase::m_dofNumber;
    using KernelBase::m_matrix;
    using KernelBase::m_rhs;
    using KernelBase::elemsToNodes;
    using KernelBase::constitutiveUpdate;
    using KernelBase::m_finiteElementSpace;

    using StackVars = StackVariables< numNodesPerElem,
                                      numNodesPerElem >;



    Kernels( arrayView1d< globalIndex const > const & inputDofNumber,
             ParallelMatrix & inputMatrix,
             ParallelVector & inputRhs,
             NodeManager const & nodeManager,
             SUBREGION_TYPE const & elementSubRegion,
             FiniteElementBase const * const finiteElementSpace,
             CONSTITUTIVE_TYPE * const inputConstitutiveType,
             Parameters const & GEOSX_UNUSED_PARAM( parameters ) )://,
      KernelBase( inputDofNumber,
                  inputMatrix,
                  inputRhs,
                  nodeManager,
                  elementSubRegion,
                  finiteElementSpace,
                  inputConstitutiveType,
                  inputConstitutiveType->createKernelWrapper() ),
      m_disp( nodeManager.totalDisplacement()),
      m_uhat( nodeManager.incrementalDisplacement()),
      dNdX( elementSubRegion.template getReference< array3d< R1Tensor > >( dataRepository::keys::dNdX )),
      detJ( elementSubRegion.template getReference< array2d< real64 > >( dataRepository::keys::detJ ) )//,
    {}

    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;
    arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

    arrayView3d< R1Tensor const > const dNdX;
    arrayView2d< real64 const > const detJ;


    template< typename STACK_VARIABLE_TYPE >
    //    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void preKernel( localIndex const k,
                    STACK_VARIABLE_TYPE & stack ) const
    {
      for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const localNodeIndex = elemsToNodes( k, a );

        stack.u_local[ a ] = m_disp[ localNodeIndex ];
        stack.uhat_local[ a ] = m_uhat[ localNodeIndex ];

        for( int i=0; i<3; ++i )
        {
          stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
        }
      }

    }

    template< typename STACK_VARIABLE_TYPE >
    //    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void updateKernel( localIndex const k,
                       localIndex const q,
                       STACK_VARIABLE_TYPE & stack ) const
    {
      real64 strainInc[6] = {0};
      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        strainInc[0] = strainInc[0] + dNdX( k, q, a )[0] * stack.uhat_local[a][0];
        strainInc[1] = strainInc[1] + dNdX( k, q, a )[1] * stack.uhat_local[a][1];
        strainInc[2] = strainInc[2] + dNdX( k, q, a )[2] * stack.uhat_local[a][2];
        strainInc[3] = strainInc[3] + dNdX( k, q, a )[2] * stack.uhat_local[a][1] +
                       dNdX( k, q, a )[1] * stack.uhat_local[a][2];
        strainInc[4] = strainInc[4] + dNdX( k, q, a )[2] * stack.uhat_local[a][0] +
                       dNdX( k, q, a )[0] * stack.uhat_local[a][2];
        strainInc[5] = strainInc[5] + dNdX( k, q, a )[1] * stack.uhat_local[a][0] +
                       dNdX( k, q, a )[0] * stack.uhat_local[a][1];
      }

      constitutiveUpdate.SmallStrain( k, q, strainInc );

      GEOSX_UNUSED_VAR( q )
      constitutiveUpdate.GetStiffness( k, stack.constitutiveStiffness );
    }


    template< typename PARAMETERS_TYPE,
              typename STACK_VARIABLE_TYPE,
              typename DYNAMICS_LAMBDA = std::function< void( localIndex, localIndex) > >
    //    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void stiffnessKernel( localIndex const k,
                          localIndex const q,
                          PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                          STACK_VARIABLE_TYPE & stack,
                          DYNAMICS_LAMBDA && dynamicsTerms = [] ( localIndex, localIndex){} ) const
    {
      for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
      {
        for( localIndex b=0; b<NUM_NODES_PER_ELEM; ++b )
        {
          real64 const (&c)[6][6] = stack.constitutiveStiffness;
          stack.localJacobian[ a*3+0 ][ b*3+0 ] -= ( c[0][0]*dNdX( k, q, a )[0]*dNdX( k, q, b )[0] +
                                                     c[5][5]*dNdX( k, q, a )[1]*dNdX( k, q, b )[1] +
                                                     c[4][4]*dNdX( k, q, a )[2]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+0 ][ b*3+1 ] -= ( c[5][5]*dNdX( k, q, a )[1]*dNdX( k, q, b )[0] +
                                                     c[0][1]*dNdX( k, q, a )[0]*dNdX( k, q, b )[1] ) * detJ( k, q );

          stack.localJacobian[ a*3+0 ][ b*3+2 ] -= ( c[4][4]*dNdX( k, q, a )[2]*dNdX( k, q, b )[0] +
                                                     c[0][2]*dNdX( k, q, a )[0]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+1 ][ b*3+1 ] -= ( c[5][5]*dNdX( k, q, a )[0]*dNdX( k, q, b )[0] +
                                                     c[1][1]*dNdX( k, q, a )[1]*dNdX( k, q, b )[1] +
                                                     c[3][3]*dNdX( k, q, a )[2]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+1 ][ b*3+0 ] -= ( c[0][1]*dNdX( k, q, a )[1]*dNdX( k, q, b )[0] +
                                                     c[5][5]*dNdX( k, q, a )[0]*dNdX( k, q, b )[1] ) * detJ( k, q );

          stack.localJacobian[ a*3+1 ][ b*3+2 ] -= ( c[3][3]*dNdX( k, q, a )[2]*dNdX( k, q, b )[1] +
                                                     c[1][2]*dNdX( k, q, a )[1]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+2 ][ b*3+0 ] -= ( c[0][2]*dNdX( k, q, a )[2]*dNdX( k, q, b )[0] +
                                                     c[4][4]*dNdX( k, q, a )[0]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+2 ][ b*3+1 ] -= ( c[1][2]*dNdX( k, q, a )[2]*dNdX( k, q, b )[1] +
                                                     c[3][3]*dNdX( k, q, a )[1]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+2 ][ b*3+2 ] -= ( c[4][4]*dNdX( k, q, a )[0]*dNdX( k, q, b )[0] +
                                                     c[3][3]*dNdX( k, q, a )[1]*dNdX( k, q, b )[1] +
                                                     c[2][2]*dNdX( k, q, a )[2]*dNdX( k, q, b )[2] ) * detJ( k, q );

          dynamicsTerms( a, b );
        }
      }
    }

    template< typename PARAMETERS_TYPE,
              typename STACK_VARIABLE_TYPE,
              typename DYNAMICS_LAMBDA = std::function< void( real64 * ) > >
    //    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void integrationKernel( localIndex const k,
                            localIndex const q,
                            PARAMETERS_TYPE const & parameters,
                            STACK_VARIABLE_TYPE & stack,
                            DYNAMICS_LAMBDA && stressModifier = [] ( real64 * ) {} ) const
    {
      real64 stress[6] = { constitutiveUpdate.m_stress( k, q, 0 ),
                           constitutiveUpdate.m_stress( k, q, 1 ),
                           constitutiveUpdate.m_stress( k, q, 2 ),
                           constitutiveUpdate.m_stress( k, q, 3 ),
                           constitutiveUpdate.m_stress( k, q, 4 ),
                           constitutiveUpdate.m_stress( k, q, 5 ) };

      stressModifier( stress );

      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        stack.localResidual[ a * 3 + 0 ] -= ( stress[ 0 ] * dNdX( k, q, a )[ 0 ] +
                                              stress[ 5 ] * dNdX( k, q, a )[ 1 ] +
                                              stress[ 4 ] * dNdX( k, q, a )[ 2 ] -
                                              parameters.m_gravityVector[0] ) * detJ( k, q );
        stack.localResidual[ a * 3 + 1 ] -= ( stress[ 5 ] * dNdX( k, q, a )[ 0 ] +
                                              stress[ 1 ] * dNdX( k, q, a )[ 1 ] +
                                              stress[ 3 ] * dNdX( k, q, a )[ 2 ] -
                                              parameters.m_gravityVector[1] ) * detJ( k, q );
        stack.localResidual[ a * 3 + 2 ] -= ( stress[ 4 ] * dNdX( k, q, a )[ 0 ] +
                                              stress[ 3 ] * dNdX( k, q, a )[ 1 ] +
                                              stress[ 2 ] * dNdX( k, q, a )[ 2 ] -
                                              parameters.m_gravityVector[2] ) * detJ( k, q );
      }
    }

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
    //    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    real64 postKernel( PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                       STACK_VARIABLE_TYPE & stack ) const
    {
      real64 meanForce = 0;
      for( localIndex a=0; a<stack.numRows; ++a )
      {
        meanForce = std::max( meanForce, stack.localResidual[a] );
        //        meanForce += fabs( stack.localResidual[a] );
      }
      //      meanForce /= stack.ndof;


      m_matrix.add( stack.localRowDofIndex,
                    stack.localColDofIndex,
                    &(stack.localJacobian[0][0]),
                    stack.numRows,
                    stack.numCols );

      m_rhs.add( stack.localRowDofIndex,
                 stack.localResidual,
                 stack.numRows );

      return meanForce;
    }

  };

};

//
//
//class ImplicitNewmark
//{
//public:
//  using Base = QuasiStatic;
//
//  struct Parameters : public Base::Parameters
//  {
//    Parameters( real64 const inputGravityVector[3],
//                real64 const inputNewmarkGamma,
//                real64 const inputNewmarkBeta,
//                real64 const inputMassDamping,
//                real64 const inputStiffnessDamping,
//                real64 const inputDt ):
//    Base::Parameters( inputGravityVector ),
//      newmarkGamma(inputNewmarkGamma),
//      newmarkBeta(inputNewmarkBeta),
//      massDamping(inputMassDamping),
//      stiffnessDamping(inputStiffnessDamping),
//      dt(inputDt)
//    {}
//
//    real64 const newmarkGamma;
//    real64 const newmarkBeta;
//    real64 const massDamping;
//    real64 const stiffnessDamping;
//    real64 const dt;
//  };
//
//  template< int NUM_NODES_PER_ELEM, int NUM_DOF_PER_NODE >
//  struct StackVariables : Base::StackVariables< NUM_NODES_PER_ELEM, NUM_DOF_PER_NODE >
//  {
//public:
//    using StackVariablesBase = Base::StackVariables< NUM_NODES_PER_ELEM, NUM_DOF_PER_NODE >;
//
//    using StackVariablesBase::numNodesPerElem;
//    using StackVariablesBase::numDofPerNode;
//    using StackVariablesBase::ndof;
//
////      GEOSX_HOST_DEVICE
//    StackVariables():
//      StackVariablesBase(),
//      dRdU_InertiaMassDamping{ {0.0} },
//      vtilde_local{ { 0.0, 0.0, 0.0} },
//      uhattilde_local{ { 0.0, 0.0, 0.0} }
//    {}
//
//    real64 dRdU_InertiaMassDamping[ ndof ][ ndof ];
//
//    R1Tensor vtilde_local[NUM_NODES_PER_ELEM];
//    R1Tensor uhattilde_local[NUM_NODES_PER_ELEM];
//  };
//
//  template< typename SUBREGION_TYPE,
//            typename CONSTITUTIVE_TYPE >
//  class Kernels : public Base::Kernels< SUBREGION_TYPE, CONSTITUTIVE_TYPE >
//  {
//  public:
//
//    using KernelBase = Base::Kernels<SUBREGION_TYPE,CONSTITUTIVE_TYPE>;
//    using KernelBase::m_dofNumber;
//    using KernelBase::m_matrix;
//    using KernelBase::m_rhs;
//    using KernelBase::elemsToNodes;
//    using KernelBase::constitutiveUpdate;
//    using KernelBase::m_disp;
//    using KernelBase::m_uhat;
//    using KernelBase::dNdX;
//    using KernelBase::detJ;
//    using KernelBase::m_finiteElementSpace;
//
//
//
//    Kernels( arrayView1d< globalIndex const > const & inputDofNumber,
//             ParallelMatrix & inputMatrix,
//             ParallelVector & inputRhs,
//             NodeManager const & nodeManager,
//             SUBREGION_TYPE const & elementSubRegion,
//             FiniteElementBase const * const finiteElementSpace,
//             CONSTITUTIVE_TYPE & constitutiveModel,
//             Base::Parameters const & parameters ):
//      KernelBase( inputDofNumber,
//                 inputMatrix,
//                 inputRhs,
//                 nodeManager,
//                 elementSubRegion,
//                 finiteElementSpace,
//                 constitutiveModel,
//                 parameters ),
//      m_vtilde(nodeManager.totalDisplacement()),
//      m_uhattilde(nodeManager.totalDisplacement()),
//      m_density(constitutiveModel.getDensity())
//    {}
//
//    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_vtilde;
//    arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhattilde;
//    arrayView2d< real64 const > const m_density;
//
//    template< typename STACK_VARIABLE_TYPE >
//  //    GEOSX_HOST_DEVICE
//    GEOSX_FORCE_INLINE
//    void preKernel( localIndex const k,
//                    STACK_VARIABLE_TYPE & stack ) const
//    {
//      for( localIndex a=0; a<STACK_VARIABLE_TYPE::numNodesPerElem; ++a )
//      {
//        localIndex const localNodeIndex = elemsToNodes( k, a );
//
//        stack.u_local[ a ] = m_disp[ localNodeIndex ];
//        stack.uhat_local[ a ] = m_uhat[ localNodeIndex ];
//        stack.vtilde_local[ a ] = m_vtilde[ localNodeIndex ];
//        stack.uhattilde_local[ a ] = m_uhattilde[ localNodeIndex ];
//
//        for( int i=0; i<3; ++i )
//        {
//          stack.elementLocalDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
//        }
//      }
//
//    }
//
//    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
//  //    GEOSX_HOST_DEVICE
//    GEOSX_FORCE_INLINE
//    void stiffnessKernel( localIndex const k,
//                          localIndex const q,
//                          PARAMETERS_TYPE const & parameters,
//                          STACK_VARIABLE_TYPE & stack ) const
//    {
//
//      std::vector< double > const & N = m_finiteElementSpace->values( q );
//
////      real64 N[STACK_VARIABLE_TYPE::numNodesPerElem];
//      real64 const & massDamping = parameters.massDamping;
//      real64 const & newmarkGamma = parameters.newmarkGamma;
//      real64 const & newmarkBeta = parameters.newmarkBeta;
//      real64 const & dt = parameters.dt;
//
//      KernelBase::stiffnessKernel( k, q, parameters, stack, [&]( localIndex const a, localIndex const b )
//      {
//        real64 integrationFactor = m_density( k, q ) * N[a] * N[b] * detJ(k,q);
//        real64 temp1 = ( massDamping * newmarkGamma/( newmarkBeta * dt ) + 1.0 / ( newmarkBeta * dt * dt ) )*
// integrationFactor;
//
//        constexpr int nsdof = STACK_VARIABLE_TYPE::numDofPerNode;
//        for( int i=0; i<nsdof; ++i )
//        {
//          realT const acc = 1.0 / ( newmarkBeta * dt * dt ) * ( stack.uhat_local[b][i] - stack.uhattilde_local[b][i]
// );
//          realT const vel = stack.vtilde_local[b][i] + newmarkGamma/( newmarkBeta * dt ) *( stack.uhat_local[b][i] -
// stack.uhattilde_local[b][i] );
//
//          stack.dRdU_InertiaMassDamping[ a*nsdof+i][ b*nsdof+i ] -= temp1;
//          stack.localResidual[ a*nsdof+i ] -= ( massDamping * vel + acc ) * integrationFactor;
//        }
//      } );
//    }
//
//    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
//  //    GEOSX_HOST_DEVICE
//    GEOSX_FORCE_INLINE
//    real64 postKernel( PARAMETERS_TYPE const & parameters,
//                       STACK_VARIABLE_TYPE & stack ) const
//    {
//      constexpr int nsdof = STACK_VARIABLE_TYPE::numDofPerNode;
//      real64 const & stiffnessDamping = parameters.stiffnessDamping;
//      real64 const & newmarkGamma     = parameters.newmarkGamma;
//      real64 const & newmarkBeta      = parameters.newmarkBeta;
//      real64 const & dt               = parameters.dt;
//
//      for( localIndex a=0; a<STACK_VARIABLE_TYPE::numNodesPerElem; ++a )
//      {
//        for( localIndex b=0; b<STACK_VARIABLE_TYPE::numNodesPerElem; ++b )
//        {
//          for( int i=0; i<nsdof; ++i )
//          {
//            for( int j=0; j<nsdof; ++j )
//            {
//              stack.localResidual[ a*nsdof+i ] += stiffnessDamping * stack.localJacobian[ a*nsdof+i][ b*nsdof+j ] *
//                                               ( stack.vtilde_local[b][j] + newmarkGamma/(newmarkBeta *
// dt)*(stack.uhat_local[b][j]-stack.uhattilde_local[b][j]) );
//
//              stack.localJacobian[a*nsdof+i][b*nsdof+j] += stack.localJacobian[a][b] * (1.0 + stiffnessDamping *
// newmarkGamma / ( newmarkBeta * dt ) ) +
//                                           stack.dRdU_InertiaMassDamping[ a ][ b ] ;
//            }
//          }
//        }
//      }
//
//      for( localIndex a=0 ; a<STACK_VARIABLE_TYPE::ndof ; ++a )
//      {
//        for( localIndex b=0 ; b<STACK_VARIABLE_TYPE::ndof ; ++b )
//        {
//          stack.localJacobian[a][b] += stack.localJacobian[a][b] * (1.0 + stiffnessDamping * newmarkGamma / (
// newmarkBeta * dt ) ) +
//                                       stack.dRdU_InertiaMassDamping[ a ][ b ] ;
//        }
//      }
//
//      return KernelBase::postKernel( parameters, stack );
//    }
//  };
//
//};




//class SmallStrainFracturePenaltyContact
//{
//public:
//  using Base = physicsLoopInterface::FiniteElementRegionLoop;
//  static constexpr int maxNumNodesPerFace = 4;
//  static constexpr int numTestDofPerSP = 3;
//  static constexpr int numTrialDofPerSP = 3;
//
//
//  template< int ,
//            int  >
//  struct StackVariables : Base::StackVariables< maxNumNodesPerFace*numTestDofPerSP*2,
//                                                maxNumNodesPerFace*numTrialDofPerSP*2 >
//  {
//public:
//    using StackVariablesBase = Base::StackVariables< maxNumNodesPerFace*numTestDofPerSP*2,
//                                                     maxNumNodesPerFace*numTrialDofPerSP*2 >;
//    using StackVariablesBase::numRows;
//    using StackVariablesBase::numCols;
//
////      GEOSX_HOST_DEVICE
//    StackVariables():
//      StackVariablesBase()
//    {}
//  };
//
//  template< typename SUBREGION_TYPE,
//            typename CONSTITUTIVE_TYPE,
//            int NUM_NODES_PER_ELEM,
//            int >
//  using SparsityKernels = Base::Kernels< SUBREGION_TYPE,
//                                         CONSTITUTIVE_TYPE,
//                                         NUM_NODES_PER_ELEM,
//                                         NUM_NODES_PER_ELEM,
//                                         numTestDofPerSP,
//                                         numTrialDofPerSP >
//  ;
//
//  template< typename SUBREGION_TYPE,
//            typename CONSTITUTIVE_TYPE,
//            int NUM_NODES_PER_ELEM,
//            int >
//  class Kernels : public Base::Kernels< SUBREGION_TYPE,
//                                        CONSTITUTIVE_TYPE,
//                                        NUM_NODES_PER_ELEM,
//                                        NUM_NODES_PER_ELEM,
//                                        numTestDofPerSP,
//                                        numTrialDofPerSP >
//  {
//public:
//    using KernelBase = Base::Kernels< SUBREGION_TYPE,
//                                      CONSTITUTIVE_TYPE,
//                                      NUM_NODES_PER_ELEM,
//                                      NUM_NODES_PER_ELEM,
//                                      numTestDofPerSP,
//                                      numTrialDofPerSP >;
//
//    static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;
//
//    using KernelBase::m_dofNumber;
//    using KernelBase::m_matrix;
//    using KernelBase::m_rhs;
//    using KernelBase::elemsToNodes;
//    using KernelBase::constitutiveUpdate;
//    using KernelBase::m_finiteElementSpace;
//
//    using StackVars = StackVariables< numNodesPerElem,
//                                      numNodesPerElem >;
//
//
//
//    Kernels( arrayView1d< globalIndex const > const & inputDofNumber,
//             ParallelMatrix & inputMatrix,
//             ParallelVector & inputRhs,
//             NodeManager const & nodeManager,
//             SUBREGION_TYPE const & elementSubRegion,
//             FiniteElementBase const * const finiteElementSpace,
//             CONSTITUTIVE_TYPE & inputConstitutiveType,
//             Parameters const & GEOSX_UNUSED_PARAM( parameters ) )://,
//      KernelBase( inputDofNumber,
//                  inputMatrix,
//                  inputRhs,
//                  nodeManager,
//                  elementSubRegion,
//                  finiteElementSpace,
//                  inputConstitutiveType,
//                  inputConstitutiveType.createKernelWrapper() ),
//      m_faceNormal( faceManager->faceNormal() ),
//      m_facesToNodes( faceManager->nodeList() ),
//      area( subRegion.getElementArea() ),
//      elemsToFaces( subRegion.faceList() ),
//
//      dNdX( elementSubRegion.template getReference< array3d< R1Tensor > >( dataRepository::keys::dNdX )),
//      detJ( elementSubRegion.template getReference< array2d< real64 > >( dataRepository::keys::detJ ) )//,
//    {}
//
//    arrayView1d< R1Tensor const > const m_faceNormal;
//    ArrayOfArraysView< localIndex const > const m_facesToNodes;
//    arrayView1d< real64 > const area;
//    arrayView2d< localIndex const > const elemsToFaces;
//
//
//
//    template< typename STACK_VARIABLE_TYPE >
//    //    GEOSX_HOST_DEVICE
//    GEOSX_FORCE_INLINE
//    void preKernel( localIndex const k,
//                    STACK_VARIABLE_TYPE & stack ) const
//    {
//      R1Tensor Nbar = m_faceNormal[elemsToFaces[kfe][0]];
//      Nbar -= m_faceNormal[elemsToFaces[kfe][1]];
//      Nbar.Normalize();
//
//      localIndex const kf0 = elemsToFaces[kfe][0];
//      localIndex const kf1 = elemsToFaces[kfe][1];
//      localIndex const numNodesPerFace=m_facesToNodes.sizeOfArray( kf0 );
//      real64 const Ja = area[kfe] / numNodesPerFace;
//
//      stackArray1d< globalIndex, maxDofPerElem > rowDOF( numNodesPerFace*3*2 );
//      stackArray1d< real64, maxDofPerElem > nodeRHS( numNodesPerFace*3*2 );
//      stackArray2d< real64, maxDofPerElem *maxDofPerElem > dRdP( numNodesPerFace*3*2, numNodesPerFace*3*2 );
//
//      for( localIndex a=0; a<numNodesPerFace; ++a )
//      {
//        localIndex const node0 = facesToNodes[kf0][a];
//        localIndex const node1 = facesToNodes[kf1][ a==0 ? a : numNodesPerFace-a ];
//
//        for( int i=0; i<3; ++i )
//        {
//          rowDOF[3*a+i]                     = nodeDofNumber[node0]+i;
//          rowDOF[3*(numNodesPerFace + a)+i] = nodeDofNumber[node1]+i;
//        }
//
//        R1Tensor gap = u[node1];
//        gap -= u[node0];
//        real64 const gapNormal = Dot( gap, Nbar );
//
//
//        if( gapNormal < 0 )
//        {
//          R1Tensor penaltyForce = Nbar;
//          penaltyForce *= -contactStiffness * gapNormal * Ja;
//          for( int i=0; i<3; ++i )
//          {
//            fc[node0] -= penaltyForce;
//            fc[node1] += penaltyForce;
//            nodeRHS[3*a+i]                     -= penaltyForce[i];
//            nodeRHS[3*(numNodesPerFace + a)+i] += penaltyForce[i];
//
//            dRdP( 3*a+i, 3*a+i )                                         -= contactStiffness * Ja * Nbar[i] * Nbar[i];
//            dRdP( 3*a+i, 3*(numNodesPerFace + a)+i )                     += contactStiffness * Ja * Nbar[i] * Nbar[i];
//            dRdP( 3*(numNodesPerFace + a)+i, 3*a+i )                     += contactStiffness * Ja * Nbar[i] * Nbar[i];
//            dRdP( 3*(numNodesPerFace + a)+i, 3*(numNodesPerFace + a)+i ) -= contactStiffness * Ja * Nbar[i] * Nbar[i];
//          }
//        }
//      }
//
//      m_residual.add( rowDOF, nodeRHS );
//      m_matrix.add( rowDOF, rowDOF, dRdP );
//
//
//    }
//
//    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
//    //    GEOSX_HOST_DEVICE
//    GEOSX_FORCE_INLINE
//    real64 postKernel( PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
//                       STACK_VARIABLE_TYPE & stack ) const
//    {
//      real64 meanForce = 0;
//      for( localIndex a=0; a<stack.numRows; ++a )
//      {
//        meanForce = std::max( meanForce, stack.localResidual[a] );
//        //        meanForce += fabs( stack.localResidual[a] );
//      }
//      //      meanForce /= stack.ndof;
//
//
//      m_matrix.add( stack.localRowDofIndex,
//                    stack.localColDofIndex,
//                    &(stack.localJacobian[0][0]),
//                    stack.numRows,
//                    stack.numCols );
//
//      m_rhs.add( stack.localRowDofIndex,
//                 stack.localResidual,
//                 stack.numRows );
//
//      return meanForce;
//    }
//  };
//};

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx
