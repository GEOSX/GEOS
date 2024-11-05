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
 * @file SinglePhaseBaseKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEKERNELS_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"

namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** MobilityKernel ********************************/

struct MobilityKernel
{
  GEOS_HOST_DEVICE
  inline
  static void
  compute( real64 const & dens,
           real64 const & dDens_dPres,
           real64 const & visc,
           real64 const & dVisc_dPres,
           real64 & mob,
           real64 & dMob_dPres )
  {
    mob = dens / visc;
    dMob_dPres = dDens_dPres / visc - mob / visc * dVisc_dPres;
  }

  GEOS_HOST_DEVICE
  inline
  static void
  compute( real64 const & dens,
           real64 const & visc,
           real64 & mob )
  {
    mob = dens / visc;
  }

  template< typename POLICY >
  static void launch( localIndex const size,
                      arrayView2d< real64 const > const & dens,
                      arrayView2d< real64 const > const & dDens_dPres,
                      arrayView2d< real64 const > const & visc,
                      arrayView2d< real64 const > const & dVisc_dPres,
                      arrayView1d< real64 > const & mob,
                      arrayView1d< real64 > const & dMob_dPres )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      compute( dens[a][0],
               dDens_dPres[a][0],
               visc[a][0],
               dVisc_dPres[a][0],
               mob[a],
               dMob_dPres[a] );
    } );
  }

  template< typename POLICY >
  static void launch( localIndex const size,
                      arrayView2d< real64 const > const & dens,
                      arrayView2d< real64 const > const & visc,
                      arrayView1d< real64 > const & mob )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      compute( dens[a][0],
               visc[a][0],
               mob[a] );
    } );
  }
};

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @brief Internal struct to provide no-op defaults used in the inclusion
 *   of lambda functions into kernel component functions.
 * @struct NoOpFunc
 */
struct NoOpFunc
{
  template< typename ... Ts >
  GEOS_HOST_DEVICE
  constexpr void
  operator()( Ts && ... ) const {}
};

/**
 * @class ElementBasedAssemblyKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation
 */
template< typename SUBREGION_TYPE, integer NUM_DOF >
class ElementBasedAssemblyKernel
{

public:

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_DOF;

  /**
   * @brief Constructor
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( globalIndex const rankOffset,
                              string const dofKey,
                              SUBREGION_TYPE const & subRegion,
                              constitutive::SingleFluidBase const & fluid,
                              constitutive::CoupledSolidBase const & solid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    :
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.template getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_volume( subRegion.getElementVolume() ),
    m_deltaVolume( subRegion.template getField< fields::flow::deltaVolume >() ),
    m_porosity( solid.getPorosity() ),
    m_dPoro_dPres( solid.getDporosity_dPressure() ),
    m_density( fluid.density() ),
    m_dDensity_dPres( fluid.dDensity_dPressure() ),
    m_mass_n( subRegion.template getField< fields::flow::mass_n >() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    // Pore volume information

    /// Pore volume at time n+1
    real64 poreVolume = 0.0;

    /// Derivative of pore volume with respect to pressure
    real64 dPoreVolume_dPres = 0.0;

    // Residual information

    /// Index of the local row corresponding to this element
    localIndex localRow = -1;

    /// Index of the matrix row/column corresponding to the dof in this element
    globalIndex dofIndices[numDof]{};

    /// Storage for the element local residual vector
    real64 localResidual[numEqn]{};

    /// Storage for the element local Jacobian matrix
    real64 localJacobian[numEqn][numDof]{};

  };

  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOS_HOST_DEVICE
  integer elemGhostRank( localIndex const ei ) const
  { return m_elemGhostRank( ei ); }


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    // initialize the pore volume
    stack.poreVolume = ( m_volume[ei] + m_deltaVolume[ei] ) * m_porosity[ei][0];
    stack.dPoreVolume_dPres = ( m_volume[ei] + m_deltaVolume[ei] ) * m_dPoro_dPres[ei][0];

    // set row index and degrees of freedom indices for this element
    stack.localRow = m_dofNumber[ei] - m_rankOffset;
    for( integer idof = 0; idof < numDof; ++idof )
    {
      stack.dofIndices[idof] = m_dofNumber[ei] + idof;
    }
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack,
                            FUNC && kernelOp = NoOpFunc{} ) const
  {
    // Residual contribution is mass conservation in the cell
    stack.localResidual[0] = stack.poreVolume * m_density[ei][0] - m_mass_n[ei];

    // Derivative of residual wrt to pressure in the cell
    stack.localJacobian[0][0] = stack.dPoreVolume_dPres * m_density[ei][0] + m_dDensity_dPres[ei][0] * stack.poreVolume;

    // Customize the kernel with this lambda
    kernelOp();
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const GEOS_UNUSED_PARAM( ei ),
                 StackVariables & stack ) const
  {
    // add contribution to global residual and jacobian (no need for atomics here)
    m_localMatrix.template addToRow< serialAtomic >( stack.localRow,
                                                     stack.dofIndices,
                                                     stack.localJacobian[0],
                                                     numDof );
    m_localRhs[stack.localRow] += stack.localResidual[0];

  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.elemGhostRank( ei ) >= 0 )
      {
        return;
      }

      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( ei, stack );
      kernelComponent.computeAccumulation( ei, stack );
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;

  /// View on the element volumes
  arrayView1d< real64 const > const m_volume;
  arrayView1d< real64 const > const m_deltaVolume;

  /// Views on the porosity
  arrayView2d< real64 const > const m_porosity;
  arrayView2d< real64 const > const m_dPoro_dPres;

  /// Views on density
  arrayView2d< real64 const > const m_density;
  arrayView2d< real64 const > const m_dDensity_dPres;

  /// View on mass
  arrayView1d< real64 const > const m_mass_n;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

};

/**
 * @class SurfaceElementBasedAssemblyKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation in SurfaceElementSubRegion
 */
class SurfaceElementBasedAssemblyKernel : public ElementBasedAssemblyKernel< SurfaceElementSubRegion, 1 >
{

public:

  using Base = ElementBasedAssemblyKernel< SurfaceElementSubRegion, 1 >;

  /**
   * @brief Constructor
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  SurfaceElementBasedAssemblyKernel( globalIndex const rankOffset,
                                     string const dofKey,
                                     SurfaceElementSubRegion const & subRegion,
                                     constitutive::SingleFluidBase const & fluid,
                                     constitutive::CoupledSolidBase const & solid,
                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                     arrayView1d< real64 > const & localRhs )
    : Base( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs )
    , m_creationMass( subRegion.getField< fields::flow::massCreated >() )
  {}

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            Base::StackVariables & stack ) const
  {
    Base::computeAccumulation( ei, stack, [&] ()
    {
      if( Base::m_mass_n[ei] > 1.1 * m_creationMass[ei] )
      {
        stack.localResidual[0] += m_creationMass[ei] * 0.25;
      }
    } );
  }

protected:

  arrayView1d< real64 const > const m_creationMass;

};

/**
 * @class ElementBasedAssemblyKernelFactory
 */
class ElementBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const dofKey,
                   CellElementSubRegion const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_DOF = 1;

    ElementBasedAssemblyKernel< CellElementSubRegion, NUM_DOF >
    kernel( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
    ElementBasedAssemblyKernel< CellElementSubRegion, NUM_DOF >::template launch< POLICY >( subRegion.size(), kernel );
  }

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const dofKey,
                   SurfaceElementSubRegion const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    SurfaceElementBasedAssemblyKernel
      kernel( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
    SurfaceElementBasedAssemblyKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< 1 >
{
public:

  using Base = physicsSolverBaseKernels::ResidualNormKernelBase< 1 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      ElementSubRegionBase const & subRegion,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_mass_n( subRegion.template getField< fields::flow::mass_n >() )
  {}

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    real64 const massNormalizer = LvArray::math::max( m_minNormalizer, m_mass_n[ei] );
    real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow] ) / massNormalizer;
    if( valMass > stack.localValue[0] )
    {
      stack.localValue[0] = valMass;
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    real64 const massNormalizer = LvArray::math::max( m_minNormalizer, m_mass_n[ei] );
    stack.localValue[0] += m_localResidual[stack.localRow] * m_localResidual[stack.localRow];
    stack.localNormalizer[0] += massNormalizer;
  }


protected:

  /// View on mass at the previous converged time step
  arrayView1d< real64 const > const m_mass_n;

};

/**
 * @class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] normType the type of norm used (Linf or L2)
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( physicsSolverBaseKernels::NormType const normType,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[1],
                   real64 (& residualNormalizer)[1] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank, subRegion, minNormalizer );
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      ResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      ResidualNormKernel::launchL2< POLICY >( subRegion.size(), kernel, residualNorm, residualNormalizer );
    }
  }

};

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY >
  static std::pair< integer, real64 > launch( arrayView1d< real64 const > const & localSolution,
                                              globalIndex const rankOffset,
                                              arrayView1d< globalIndex const > const & dofNumber,
                                              arrayView1d< integer const > const & ghostRank,
                                              arrayView1d< real64 const > const & pres,
                                              real64 const scalingFactor )
  {
    RAJA::ReduceSum< ReducePolicy< POLICY >, integer > numNegativePressures( 0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minPres( 0.0 );

    forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 && dofNumber[ei] >= 0 )
      {
        localIndex const lid = dofNumber[ei] - rankOffset;
        real64 const newPres = pres[ei] + scalingFactor * localSolution[lid];

        if( newPres < 0.0 )
        {
          numNegativePressures += 1;
          minPres.min( newPres );
        }
      }

    } );

    return { numNegativePressures.get(), minPres.get() };
  }

};

/******************************** ScalingForSystemSolutionKernel ********************************/

struct ScalingForSystemSolutionKernel
{
  template< typename POLICY >
  static std::pair< real64, real64 > launch( arrayView1d< real64 const > const & localSolution,
                                             globalIndex const rankOffset,
                                             arrayView1d< globalIndex const > const & dofNumber,
                                             arrayView1d< integer const > const & ghostRank,
                                             real64 const maxAbsolutePresChange )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > scalingFactor( 1.0 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaPres( 0.0 );

    forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei ) mutable
    {
      if( ghostRank[ei] < 0 && dofNumber[ei] >= 0 )
      {
        localIndex const lid = dofNumber[ei] - rankOffset;

        // compute the change in pressure
        real64 const absPresChange = LvArray::math::abs( localSolution[lid] );
        maxDeltaPres.max( absPresChange );

        // maxAbsolutePresChange <= 0.0 means that scaling is disabled, and we are only collecting maxDeltaPres in this kernel
        if( maxAbsolutePresChange > 0.0 && absPresChange > maxAbsolutePresChange )
        {
          real64 const presScalingFactor = maxAbsolutePresChange / absPresChange;
          scalingFactor.min( presScalingFactor );
        }
      }

    } );

    return { scalingFactor.get(), maxDeltaPres.get() };
  }

};

/******************************** StatisticsKernel ********************************/

struct StatisticsKernel
{
  static void
  saveDeltaPressure( localIndex const size,
                     arrayView1d< real64 const > const & pres,
                     arrayView1d< real64 const > const & initPres,
                     arrayView1d< real64 > const & deltaPres )
  {
    forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      deltaPres[ei] = pres[ei] - initPres[ei];
    } );
  }

  static void
  launch( localIndex const size,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & deltaPres,
          arrayView1d< real64 const > const & temp,
          arrayView1d< real64 const > const & refPorosity,
          arrayView2d< real64 const > const & porosity,
          arrayView2d< real64 const > const & density,
          real64 & minPres,
          real64 & avgPresNumerator,
          real64 & maxPres,
          real64 & minDeltaPres,
          real64 & maxDeltaPres,
          real64 & minTemp,
          real64 & avgTempNumerator,
          real64 & maxTemp,
          real64 & totalUncompactedPoreVol,
          real64 & totalPoreVol,
          real64 & totalMass )
  {
    RAJA::ReduceMin< parallelDeviceReduce, real64 > subRegionMinPres( LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionAvgPresNumerator( 0.0 );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxPres( -LvArray::NumericLimits< real64 >::max );

    RAJA::ReduceMin< parallelDeviceReduce, real64 > subRegionMinDeltaPres( LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxDeltaPres( -LvArray::NumericLimits< real64 >::max );

    RAJA::ReduceMin< parallelDeviceReduce, real64 > subRegionMinTemp( LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionAvgTempNumerator( 0.0 );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxTemp( -LvArray::NumericLimits< real64 >::max );

    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionTotalUncompactedPoreVol( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionTotalPoreVol( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionTotalMass( 0.0 );

    forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( elemGhostRank[ei] >= 0 )
      {
        return;
      }

      // To match our "reference", we have to use reference porosity here, not the actual porosity when we compute averages
      real64 const uncompactedPoreVol = volume[ei] * refPorosity[ei];
      real64 const dynamicPoreVol = volume[ei] * porosity[ei][0];

      subRegionMinPres.min( pres[ei] );
      subRegionAvgPresNumerator += uncompactedPoreVol * pres[ei];
      subRegionMaxPres.max( pres[ei] );

      subRegionMinDeltaPres.min( deltaPres[ei] );
      subRegionMaxDeltaPres.max( deltaPres[ei] );

      subRegionMinTemp.min( temp[ei] );
      subRegionAvgTempNumerator += uncompactedPoreVol * temp[ei];
      subRegionMaxTemp.max( temp[ei] );

      subRegionTotalUncompactedPoreVol += uncompactedPoreVol;
      subRegionTotalPoreVol += dynamicPoreVol;
      subRegionTotalMass += dynamicPoreVol * density[ei][0];
    } );

    minPres = subRegionMinPres.get();
    avgPresNumerator = subRegionAvgPresNumerator.get();
    maxPres = subRegionMaxPres.get();

    minDeltaPres = subRegionMinDeltaPres.get();
    maxDeltaPres = subRegionMaxDeltaPres.get();

    minTemp = subRegionMinTemp.get();
    avgTempNumerator = subRegionAvgTempNumerator.get();
    maxTemp = subRegionMaxTemp.get();

    totalUncompactedPoreVol = subRegionTotalUncompactedPoreVol.get();
    totalPoreVol = subRegionTotalPoreVol.get();
    totalMass = subRegionTotalMass.get();
  }
};


/******************************** HydrostaticPressureKernel ********************************/

struct HydrostaticPressureKernel
{

  template< typename FLUID_WRAPPER >
  static bool
  computeHydrostaticPressure( integer const maxNumEquilIterations,
                              real64 const & equilTolerance,
                              real64 const (&gravVector)[ 3 ],
                              FLUID_WRAPPER fluidWrapper,
                              real64 const & refElevation,
                              real64 const & refPres,
                              real64 const & refDens,
                              real64 const & newElevation,
                              real64 & newPres,
                              real64 & newDens )
  {
    // Step 1: guess the pressure with the refDensity

    real64 const gravCoef = gravVector[2] * ( refElevation - newElevation );
    real64 pres0 = refPres - refDens * gravCoef;
    real64 pres1 = 0.0;

    // Step 2: compute the mass density at this elevation using the guess, and update pressure

    real64 dens = 0.0;
    real64 visc = 0.0;
    constitutive::SingleFluidBaseUpdate::computeValues( fluidWrapper,
                                                        pres0,
                                                        dens,
                                                        visc );
    pres1 = refPres - 0.5 * ( refDens + dens ) * gravCoef;

    // Step 3: fixed-point iteration until convergence

    bool equilHasConverged = false;
    for( localIndex eqIter = 0; eqIter < maxNumEquilIterations; ++eqIter )
    {

      // check convergence
      equilHasConverged = ( LvArray::math::abs( pres0 - pres1 ) < equilTolerance );
      pres0 = pres1;

      // if converged, move on
      if( equilHasConverged )
      {
        break;
      }

      // compute the density at this elevation using the previous pressure, and compute the new pressure
      constitutive::SingleFluidBaseUpdate::computeValues( fluidWrapper,
                                                          pres0,
                                                          dens,
                                                          visc );
      pres1 = refPres - 0.5 * ( refDens + dens ) * gravCoef;
    }

    // Step 4: save the hydrostatic pressure and the corresponding density

    newPres = pres1;
    newDens = dens;

    return equilHasConverged;
  }


  template< typename FLUID_WRAPPER >
  static bool
  launch( localIndex const size,
          integer const maxNumEquilIterations,
          real64 const equilTolerance,
          real64 const (&gravVector)[ 3 ],
          real64 const & minElevation,
          real64 const & elevationIncrement,
          real64 const & datumElevation,
          real64 const & datumPres,
          FLUID_WRAPPER fluidWrapper,
          arrayView1d< arrayView1d< real64 > const > elevationValues,
          arrayView1d< real64 > pressureValues )
  {
    bool hasConverged = true;

    // Step 1: compute the mass density at the datum elevation

    real64 datumDens = 0.0;
    real64 datumVisc = 0.0;

    constitutive::SingleFluidBaseUpdate::computeValues( fluidWrapper,
                                                        datumPres,
                                                        datumDens,
                                                        datumVisc );

    // Step 2: find the closest elevation to datumElevation

    forAll< parallelHostPolicy >( size, [=] ( localIndex const i )
    {
      real64 const elevation = minElevation + i * elevationIncrement;
      elevationValues[0][i] = elevation;
    } );
    integer const iRef = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                 elevationValues[0].size(),
                                                                 datumElevation );


    // Step 3: compute the mass density and pressure at the reference elevation

    array1d< real64 > dens( pressureValues.size() );

    bool const hasConvergedRef =
      computeHydrostaticPressure( maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluidWrapper,
                                  datumElevation,
                                  datumPres,
                                  datumDens,
                                  elevationValues[0][iRef],
                                  pressureValues[iRef],
                                  dens[iRef] );
    if( !hasConvergedRef )
    {
      return false;
    }


    // Step 4: for each elevation above the reference elevation, compute the pressure

    localIndex const numEntriesAboveRef = size - iRef - 1;
    forAll< serialPolicy >( numEntriesAboveRef, [=, &hasConverged] ( localIndex const i )
    {
      bool const hasConvergedAboveRef =
        computeHydrostaticPressure( maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    elevationValues[0][iRef+i],
                                    pressureValues[iRef+i],
                                    dens[iRef+i],
                                    elevationValues[0][iRef+i+1],
                                    pressureValues[iRef+i+1],
                                    dens[iRef+i+1] );
      if( !hasConvergedAboveRef )
      {
        hasConverged = false;
      }


    } );

    // Step 5: for each elevation below the reference elevation, compute the pressure

    localIndex const numEntriesBelowRef = iRef;
    forAll< serialPolicy >( numEntriesBelowRef, [=, &hasConverged] ( localIndex const i )
    {
      bool const hasConvergedBelowRef =
        computeHydrostaticPressure( maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    elevationValues[0][iRef-i],
                                    pressureValues[iRef-i],
                                    dens[iRef-i],
                                    elevationValues[0][iRef-i-1],
                                    pressureValues[iRef-i-1],
                                    dens[iRef-i-1] );
      if( !hasConvergedBelowRef )
      {
        hasConverged = false;
      }
    } );

    return hasConverged;
  }
};


} // namespace singlePhaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEKERNELS_HPP
