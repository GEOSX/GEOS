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
 * @file IsothermalCompositionalMultiphaseBaseKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEBASEKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEBASEKERNELS_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "functions/TableFunction.hpp"
#include "mesh/ElementSubRegionBase.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/KernelLaunchSelectors.hpp"

namespace geos
{


namespace isothermalCompositionalMultiphaseBaseKernels
{

static constexpr real64 minDensForDivision = 1e-10;

enum class ElementBasedAssemblyKernelFlags
{
  SimpleAccumulation = 1 << 0, // 1
  TotalMassEquation = 1 << 1, // 2
  /// Add more flags like that if needed:
  // Flag3 = 1 << 2, // 4
  // Flag4 = 1 << 3, // 8
  // Flag5 = 1 << 4, // 16
  // Flag6 = 1 << 5, // 32
  // Flag7 = 1 << 6, // 64
  // Flag8 = 1 << 7  //128
};

/******************************** PropertyKernelBase ********************************/

/**
 * @class PropertyKernelBase
 * @tparam NUM_COMP number of fluid components
 * @brief Define the base interface for the property update kernels
 */
template< integer NUM_COMP >
class PropertyKernelBase
{
public:

  /// Compile time value for the number of components
  static constexpr integer numComp = NUM_COMP;

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      kernelComponent.compute( ei );
    } );
  }

  /**
   * @brief Performs the kernel launch on a sorted array
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] targetSet the indices of the elements in which we compute the property
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          KERNEL_TYPE const & kernelComponent )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localIndex const ei = targetSet[ i ];
      kernelComponent.compute( ei );
    } );
  }

};

namespace internal
{



template< typename T, typename LAMBDA >
void kernelLaunchSelectorCompSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorCompSwitch: type should be integral" );

  switch( value )
  {
    case 1:
    { lambda( std::integral_constant< T, 1 >() ); return; }
    case 2:
    { lambda( std::integral_constant< T, 2 >() ); return; }
    case 3:
    { lambda( std::integral_constant< T, 3 >() ); return; }
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return; }
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return; }
    default:
    { GEOS_ERROR( "Unsupported number of components: " << value ); }
  }
}

} // namespace internal


/******************************** GlobalComponentFractionKernel ********************************/

/**
 * @class GlobalComponentFractionKernel
 * @tparam NUM_COMP number of fluid components
 * @brief Define the interface for the update kernel in charge of computing the phase volume fractions
 */
template< integer NUM_COMP >
class GlobalComponentFractionKernel : public PropertyKernelBase< NUM_COMP >
{
public:

  using Base = PropertyKernelBase< NUM_COMP >;
  using Base::numComp;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  GlobalComponentFractionKernel( ObjectManagerBase & subRegion )
    : Base(),
    m_compDens( subRegion.getField< fields::flow::globalCompDensity >() ),
    m_compFrac( subRegion.getField< fields::flow::globalCompFraction >() ),
    m_dCompFrac_dCompDens( subRegion.getField< fields::flow::dGlobalCompFraction_dGlobalCompDensity >() )
  {}

  /**
   * @brief Compute the phase volume fractions in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseVolFractionKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                FUNC && compFractionKernelOp = NoOpFunc{} ) const
  {
    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const compDens = m_compDens[ei];
    arraySlice1d< real64, compflow::USD_COMP - 1 > const compFrac = m_compFrac[ei];
    arraySlice2d< real64, compflow::USD_COMP_DC - 1 > const dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];

    real64 totalDensity = 0.0;

    for( integer ic = 0; ic < numComp; ++ic )
    {
      totalDensity += compDens[ic];
    }

    real64 const totalDensityInv = 1.0 / totalDensity;

    for( integer ic = 0; ic < numComp; ++ic )
    {
      compFrac[ic] = compDens[ic] * totalDensityInv;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dCompFrac_dCompDens[ic][jc] = -compFrac[ic] * totalDensityInv;
      }
      dCompFrac_dCompDens[ic][ic] += totalDensityInv;
    }

    compFractionKernelOp( compFrac, dCompFrac_dCompDens );
  }

protected:

  // inputs

  // Views on component densities
  arrayView2d< real64 const, compflow::USD_COMP > m_compDens;

  // outputs

  // Views on component fraction
  arrayView2d< real64, compflow::USD_COMP > m_compFrac;
  arrayView3d< real64, compflow::USD_COMP_DC > m_dCompFrac_dCompDens;

};

/**
 * @class GlobalComponentFractionKernelFactory
 */
class GlobalComponentFractionKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   ObjectManagerBase & subRegion )
  {
    internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      GlobalComponentFractionKernel< NUM_COMP > kernel( subRegion );
      GlobalComponentFractionKernel< NUM_COMP >::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }

};

/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @class PhaseVolumeFractionKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase volume fractions
 */
template< integer NUM_COMP, integer NUM_PHASE >
class PhaseVolumeFractionKernel : public PropertyKernelBase< NUM_COMP >
{
public:

  using Base = PropertyKernelBase< NUM_COMP >;
  using Base::numComp;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NUM_PHASE;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  PhaseVolumeFractionKernel( ObjectManagerBase & subRegion,
                             constitutive::MultiFluidBase const & fluid )
    : Base(),
    m_phaseVolFrac( subRegion.getField< fields::flow::phaseVolumeFraction >() ),
    m_dPhaseVolFrac( subRegion.getField< fields::flow::dPhaseVolumeFraction >() ),
    m_compDens( subRegion.getField< fields::flow::globalCompDensity >() ),
    m_dCompFrac_dCompDens( subRegion.getField< fields::flow::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseFrac( fluid.phaseFraction() ),
    m_dPhaseFrac( fluid.dPhaseFraction() ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() )
  {}

  /**
   * @brief Compute the phase volume fractions in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseVolFractionKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  real64 compute( localIndex const ei,
                  FUNC && phaseVolFractionKernelOp = NoOpFunc{} ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const compDens = m_compDens[ei];
    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseDens = m_phaseDens[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseDens = m_dPhaseDens[ei][0];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseFrac = m_phaseFrac[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseFrac = m_dPhaseFrac[ei][0];
    arraySlice1d< real64, compflow::USD_PHASE - 1 > const phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dPhaseVolFrac[ei];

    real64 work[numComp]{};

    // compute total density from component partial densities
    real64 totalDensity = 0.0;
    real64 const dTotalDens_dCompDens = 1.0;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      totalDensity += compDens[ic];
    }

    real64 maxDeltaPhaseVolFrac = 0.0;

    for( integer ip = 0; ip < numPhase; ++ip )
    {

      // set the saturation to zero if the phase is absent
      bool const phaseExists = (phaseFrac[ip] > 0);
      if( !phaseExists )
      {
        phaseVolFrac[ip] = 0.;
        for( integer jc = 0; jc < numComp+2; ++jc )
        {
          dPhaseVolFrac[ip][jc] = 0.;
        }
        continue;
      }

      // Expression for volume fractions: S_p = (nu_p / rho_p) * rho_t
      real64 const phaseDensInv = 1.0 / phaseDens[ip];

      // store old saturation to compute change later
      real64 const satOld = phaseVolFrac[ip];

      // compute saturation and derivatives except multiplying by the total density
      phaseVolFrac[ip] = phaseFrac[ip] * phaseDensInv;

      dPhaseVolFrac[ip][Deriv::dP] =
        (dPhaseFrac[ip][Deriv::dP] - phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dP]) * phaseDensInv;

      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseVolFrac[ip][Deriv::dC+jc] =
          (dPhaseFrac[ip][Deriv::dC+jc] - phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dC+jc]) * phaseDensInv;
      }

      // apply chain rule to convert derivatives from global component fractions to densities
      applyChainRuleInPlace( numComp, dCompFrac_dCompDens, dPhaseVolFrac[ip], work, Deriv::dC );

      // call the lambda in the phase loop to allow the reuse of the phaseVolFrac and totalDensity
      // possible use: assemble the derivatives wrt temperature
      phaseVolFractionKernelOp( ip, phaseVolFrac[ip], phaseDensInv, totalDensity );

      // now finalize the computation by multiplying by total density
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseVolFrac[ip][Deriv::dC+jc] *= totalDensity;
        dPhaseVolFrac[ip][Deriv::dC+jc] += phaseVolFrac[ip] * dTotalDens_dCompDens;
      }

      phaseVolFrac[ip] *= totalDensity;
      dPhaseVolFrac[ip][Deriv::dP] *= totalDensity;

      real64 const deltaPhaseVolFrac = LvArray::math::abs( phaseVolFrac[ip] - satOld );
      if( maxDeltaPhaseVolFrac < deltaPhaseVolFrac )
      {
        maxDeltaPhaseVolFrac = deltaPhaseVolFrac;
      }
    }
    return maxDeltaPhaseVolFrac;
  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static real64
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaPhaseVolFrac( 0.0 );
    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      real64 const deltaPhaseVolFrac = kernelComponent.compute( ei );
      maxDeltaPhaseVolFrac.max( deltaPhaseVolFrac );
    } );
    return maxDeltaPhaseVolFrac.get();
  }

protected:

  // outputs

  /// Views on phase volume fractions
  arrayView2d< real64, compflow::USD_PHASE > m_phaseVolFrac;
  arrayView3d< real64, compflow::USD_PHASE_DC > m_dPhaseVolFrac;

  // inputs

  /// Views on component densities
  arrayView2d< real64 const, compflow::USD_COMP > m_compDens;
  arrayView3d< real64 const, compflow::USD_COMP_DC > m_dCompFrac_dCompDens;

  /// Views on phase fractions
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseFrac;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseFrac;

  /// Views on phase densities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseDens;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseDens;

};

/**
 * @class PhaseVolumeFractionKernelFactory
 */
class PhaseVolumeFractionKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  template< typename POLICY >
  static real64
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   constitutive::MultiFluidBase const & fluid )
  {
    real64 maxDeltaPhaseVolFrac = 0.0;
    if( numPhase == 2 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseVolumeFractionKernel< NUM_COMP, 2 > kernel( subRegion, fluid );
        maxDeltaPhaseVolFrac = PhaseVolumeFractionKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseVolumeFractionKernel< NUM_COMP, 3 > kernel( subRegion, fluid );
        maxDeltaPhaseVolFrac = PhaseVolumeFractionKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    return maxDeltaPhaseVolFrac;
  }
};


/******************************** RelativePermeabilityUpdateKernel ********************************/

struct RelativePermeabilityUpdateKernel
{
  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  launch( localIndex const size,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }
};

/******************************** CapillaryPressureUpdateKernel ********************************/

struct CapillaryPressureUpdateKernel
{
  template< typename POLICY, typename CAPPRES_WRAPPER >
  static void
  launch( localIndex const size,
          CAPPRES_WRAPPER const & capPresWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < capPresWrapper.numGauss(); ++q )
      {
        capPresWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename CAPPRES_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          CAPPRES_WRAPPER const & capPresWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < capPresWrapper.numGauss(); ++q )
      {
        capPresWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }
};

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of accumulation and volume balance
 */
template< integer NUM_COMP, integer NUM_DOF >
class ElementBasedAssemblyKernel
{
public:

  /// Compile time value for the number of components
  static constexpr integer numComp = NUM_COMP;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_DOF;

  /**
   * @brief Constructor
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( localIndex const numPhases,
                              globalIndex const rankOffset,
                              string const dofKey,
                              ElementSubRegionBase const & subRegion,
                              constitutive::MultiFluidBase const & fluid,
                              constitutive::CoupledSolidBase const & solid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs,
                              BitFlags< ElementBasedAssemblyKernelFlags > const kernelFlags )
    : m_numPhases( numPhases ),
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_volume( subRegion.getElementVolume() ),
    m_porosity( solid.getPorosity() ),
    m_dPoro_dPres( solid.getDporosity_dPressure() ),
    m_dCompFrac_dCompDens( subRegion.getField< fields::flow::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseVolFrac( subRegion.getField< fields::flow::phaseVolumeFraction >() ),
    m_dPhaseVolFrac( subRegion.getField< fields::flow::dPhaseVolumeFraction >() ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() ),
    m_phaseCompFrac( fluid.phaseCompFraction() ),
    m_dPhaseCompFrac( fluid.dPhaseCompFraction() ),
    m_compDens( subRegion.getField< fields::flow::globalCompDensity >() ),
    m_compAmount_n( subRegion.getField< fields::flow::compAmount_n >() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs ),
    m_kernelFlags( kernelFlags )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    // Pore volume information (used by both accumulation and volume balance)

    /// Pore volume at time n+1
    real64 poreVolume = 0.0;

    /// Derivative of pore volume with respect to pressure
    real64 dPoreVolume_dPres = 0.0;

    // Residual information

    /// Index of the local row corresponding to this element
    localIndex localRow = -1;

    /// Indices of the matrix rows/columns corresponding to the dofs in this element
    globalIndex dofIndices[numDof]{};

    /// C-array storage for the element local residual vector (all equations except volume balance)
    real64 localResidual[numEqn]{};

    /// C-array storage for the element local Jacobian matrix (all equations except volume balance, all dofs)
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
    stack.poreVolume = m_volume[ei] * m_porosity[ei][0];
    stack.dPoreVolume_dPres = m_volume[ei] * m_dPoro_dPres[ei][0];

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
   * @param[in] phaseAmountKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack,
                            FUNC && phaseAmountKernelOp = NoOpFunc{} ) const
  {
    if( m_kernelFlags.isSet( ElementBasedAssemblyKernelFlags::SimpleAccumulation ) )
    {
      // ic - index of component whose conservation equation is assembled
      // (i.e. row number in local matrix)
      for( integer ic = 0; ic < numComp; ++ic )
      {
        real64 const compAmount = stack.poreVolume * m_compDens[ei][ic];
        real64 const compAmount_n = m_compAmount_n[ei][ic];

        stack.localResidual[ic] += compAmount - compAmount_n;

        // Pavel: commented below is some experiment, needs to be re-tested
        //real64 const compDens = (ic == 0 && m_compDens[ei][ic] < 1e-6) ? 1e-3 : m_compDens[ei][ic];
        real64 const dCompAmount_dP = stack.dPoreVolume_dPres * m_compDens[ei][ic];
        stack.localJacobian[ic][0] += dCompAmount_dP;

        real64 const dCompAmount_dC = stack.poreVolume;
        stack.localJacobian[ic][ic + 1] += dCompAmount_dC;
      }
    }
    else
    {
      using Deriv = constitutive::multifluid::DerivativeOffset;

      // construct the slices for variables accessed multiple times
      arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];

      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
      arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];

      arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > phaseDens = m_phaseDens[ei][0];
      arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > dPhaseDens = m_dPhaseDens[ei][0];

      arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > phaseCompFrac = m_phaseCompFrac[ei][0];
      arraySlice3d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac = m_dPhaseCompFrac[ei][0];

      // temporary work arrays
      real64 dPhaseAmount_dC[numComp]{};
      real64 dPhaseCompFrac_dC[numComp]{};

      // start with old time step values
      for( integer ic = 0; ic < numComp; ++ic )
      {
        stack.localResidual[ic] = -m_compAmount_n[ei][ic];
      }

      // sum contributions to component accumulation from each phase
      for( integer ip = 0; ip < m_numPhases; ++ip )
      {
        real64 const phaseAmount = stack.poreVolume * phaseVolFrac[ip] * phaseDens[ip];

        real64 const dPhaseAmount_dP = stack.dPoreVolume_dPres * phaseVolFrac[ip] * phaseDens[ip]
                                       + stack.poreVolume * ( dPhaseVolFrac[ip][Deriv::dP] * phaseDens[ip]
                                                              + phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dP] );

        // assemble density dependence
        applyChainRule( numComp, dCompFrac_dCompDens, dPhaseDens[ip], dPhaseAmount_dC, Deriv::dC );
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                                + phaseDens[ip] * dPhaseVolFrac[ip][Deriv::dC + jc];
          dPhaseAmount_dC[jc] *= stack.poreVolume;
        }

        // ic - index of component whose conservation equation is assembled
        // (i.e. row number in local matrix)
        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 const phaseCompAmount = phaseAmount * phaseCompFrac[ip][ic];

          real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ip][ic]
                                             + phaseAmount * dPhaseCompFrac[ip][ic][Deriv::dP];

          stack.localResidual[ic] += phaseCompAmount;
          stack.localJacobian[ic][0] += dPhaseCompAmount_dP;

          // jc - index of component w.r.t. whose compositional var the derivative is being taken
          // (i.e. col number in local matrix)

          // assemble phase composition dependence
          applyChainRule( numComp, dCompFrac_dCompDens, dPhaseCompFrac[ip][ic], dPhaseCompFrac_dC, Deriv::dC );
          for( integer jc = 0; jc < numComp; ++jc )
          {
            real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmount
                                               + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];

            stack.localJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
          }
        }

        // call the lambda in the phase loop to allow the reuse of the phase amounts and their derivatives
        // possible use: assemble the derivatives wrt temperature, and the accumulation term of the energy equation for this phase
        phaseAmountKernelOp( ip, phaseAmount, dPhaseAmount_dP, dPhaseAmount_dC );

      }
    }
  }

  /**
   * @brief Compute the local volume balance contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] phaseVolFractionSumKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeVolumeBalance( localIndex const ei,
                             StackVariables & stack,
                             FUNC && phaseVolFractionSumKernelOp = NoOpFunc{} ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];

    real64 oneMinusPhaseVolFracSum = 1.0;

    // sum contributions to component accumulation from each phase
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      oneMinusPhaseVolFracSum -= phaseVolFrac[ip];
      stack.localJacobian[numComp][0] -= dPhaseVolFrac[ip][Deriv::dP];

      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.localJacobian[numComp][jc+1] -= dPhaseVolFrac[ip][Deriv::dC+jc];
      }
    }

    // call the lambda in the phase loop to allow the reuse of the phase amounts and their derivatives
    // possible use: assemble the derivatives wrt temperature, and use oneMinusPhaseVolFracSum if poreVolume depends on temperature
    phaseVolFractionSumKernelOp( oneMinusPhaseVolFracSum );

    // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
    stack.localResidual[numComp] = stack.poreVolume * oneMinusPhaseVolFracSum;
    for( integer idof = 0; idof < numDof; ++idof )
    {
      stack.localJacobian[numComp][idof] *= stack.poreVolume;
    }
    stack.localJacobian[numComp][0] += stack.dPoreVolume_dPres * oneMinusPhaseVolFracSum;
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
    using namespace compositionalMultiphaseUtilities;

    if( m_kernelFlags.isSet( ElementBasedAssemblyKernelFlags::TotalMassEquation ) )
    {
      // apply equation/variable change transformation to the component mass balance equations
      real64 work[numDof]{};
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numDof, stack.localJacobian, work );
      shiftElementsAheadByOneAndReplaceFirstElementWithSum( numComp, stack.localResidual );
    }

    // add contribution to residual and jacobian into:
    // - the component mass balance equations (i = 0 to i = numComp-1)
    // - the volume balance equations (i = numComp)
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    integer const numRows = numComp+1;
    for( integer i = 0; i < numRows; ++i )
    {
      m_localRhs[stack.localRow + i] += stack.localResidual[i];
      m_localMatrix.addToRow< serialAtomic >( stack.localRow + i,
                                              stack.dofIndices,
                                              stack.localJacobian[i],
                                              numDof );
    }
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
      kernelComponent.computeVolumeBalance( ei, stack );
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

  /// Number of fluid phases
  integer const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;

  /// View on the element volumes
  arrayView1d< real64 const > const m_volume;

  /// Views on the porosity
  arrayView2d< real64 const > const m_porosity;
  arrayView2d< real64 const > const m_dPoro_dPres;

  /// Views on the derivatives of comp fractions wrt component density
  arrayView3d< real64 const, compflow::USD_COMP_DC > const m_dCompFrac_dCompDens;

  /// Views on the phase volume fractions
  arrayView2d< real64 const, compflow::USD_PHASE > const m_phaseVolFrac;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > const m_dPhaseVolFrac;

  /// Views on the phase densities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseDens;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const m_dPhaseDens;

  /// Views on the phase component fraction
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > const m_phaseCompFrac;
  arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > const m_dPhaseCompFrac;

  // View on component densities
  arrayView2d< real64 const, compflow::USD_COMP > m_compDens;

  // View on component amount (mass/moles) from previous time step
  arrayView2d< real64 const, compflow::USD_COMP > m_compAmount_n;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  BitFlags< ElementBasedAssemblyKernelFlags > const m_kernelFlags;
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
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
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
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   integer const useTotalMassEquation,
                   integer const useSimpleAccumulation,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   constitutive::MultiFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC()+1;

      BitFlags< ElementBasedAssemblyKernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( ElementBasedAssemblyKernelFlags::TotalMassEquation );
      if( useSimpleAccumulation )
        kernelFlags.set( ElementBasedAssemblyKernelFlags::SimpleAccumulation );

      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
      kernel( numPhases, rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs, kernelFlags );
      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }

};

/******************************** ScalingForSystemSolutionKernel ********************************/

/**
 * @class ScalingAndCheckingSystemSolutionKernelBase
 * @brief Define the kernel for scaling the solution and check its validity
 */
template< typename TYPE >
class ScalingAndCheckingSystemSolutionKernelBase
{
public:

  /**
   * @brief Create a new kernel instance
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @param[in] pressure the pressure vector
   * @param[in] compDens the component density vector
   * @param[in] pressureScalingFactor the pressure local scaling factor
   * @param[in] compDensScalingFactor the component local scaling factor
   */
  ScalingAndCheckingSystemSolutionKernelBase( globalIndex const rankOffset,
                                              integer const numComp,
                                              string const dofKey,
                                              ElementSubRegionBase const & subRegion,
                                              arrayView1d< real64 const > const localSolution,
                                              arrayView1d< real64 const > const pressure,
                                              arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                                              arrayView1d< real64 > pressureScalingFactor,
                                              arrayView1d< real64 > compDensScalingFactor )
    : m_rankOffset( rankOffset ),
    m_numComp( numComp ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_ghostRank( subRegion.ghostRank() ),
    m_localSolution( localSolution ),
    m_pressure( pressure ),   // not passed with fields::flow to be able to reuse this for wells
    m_compDens( compDens ),   // same here
    m_pressureScalingFactor( pressureScalingFactor ),
    m_compDensScalingFactor( compDensScalingFactor )
  { }

  /**
   * @struct StackVariables
   * @brief Kernel variables located on the stack
   */
  struct StackVariables
  {
    GEOS_HOST_DEVICE
    StackVariables()
    { }

    StackVariables( real64 _localMinVal )
      :
      localMinVal( _localMinVal )
    { }

    /// Index of the local row corresponding to this element
    localIndex localRow;

    /// The local value
    TYPE localMinVal;
  };

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    stack.localMinVal = 1;

    // set row index and degrees of freedom indices for this element
    stack.localRow = m_dofNumber[ei] - m_rankOffset;
  }

  /**
   * @brief Getter for the ghost rank
   * @param[in] i the looping index of the element/node/face
   * @return the ghost rank of the element/node/face
   */
  GEOS_HOST_DEVICE
  integer ghostRank( localIndex const i ) const
  { return m_ghostRank( i ); }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static TYPE
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, TYPE > minVal( 1 );
    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.ghostRank( ei ) >= 0 )
      {
        return;
      }

      StackVariables stack;
      kernelComponent.setup( ei, stack );
      kernelComponent.compute( ei, stack );
      minVal.min( stack.localMinVal );
    } );

    return minVal.get();
  }

protected:

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Number of components
  real64 const m_numComp;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_ghostRank;

  /// View on the local residual
  arrayView1d< real64 const > const m_localSolution;

  /// View on the primary variables
  arrayView1d< real64 const > const m_pressure;
  arrayView2d< real64 const, compflow::USD_COMP > const m_compDens;

  /// View on the scaling factors
  arrayView1d< real64 > const m_pressureScalingFactor;
  arrayView1d< real64 > const m_compDensScalingFactor;

};

/**
 * @class ScalingForSystemSolutionKernel
 * @brief Define the kernel for scaling the Newton update
 */
class ScalingForSystemSolutionKernel : public ScalingAndCheckingSystemSolutionKernelBase< real64 >
{
public:

  using Base = ScalingAndCheckingSystemSolutionKernelBase< real64 >;
  using Base::m_rankOffset;
  using Base::m_numComp;
  using Base::m_dofNumber;
  using Base::m_ghostRank;
  using Base::m_localSolution;
  using Base::m_pressure;
  using Base::m_compDens;
  using Base::m_pressureScalingFactor;
  using Base::m_compDensScalingFactor;

  /**
   * @brief Create a new kernel instance
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxAbsolutePresChange the max allowed absolute pressure change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] maxRelativeCompDensChange the max allowed comp density change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @param[in] pressure the pressure vector
   * @param[in] compDens the component density vector
   * @param[in] pressureScalingFactor the pressure local scaling factor
   * @param[in] compDensScalingFactor the component density local scaling factor
   */
  ScalingForSystemSolutionKernel( real64 const maxRelativePresChange,
                                  real64 const maxAbsolutePresChange,
                                  real64 const maxCompFracChange,
                                  real64 const maxRelativeCompDensChange,
                                  globalIndex const rankOffset,
                                  integer const numComp,
                                  string const dofKey,
                                  ElementSubRegionBase const & subRegion,
                                  arrayView1d< real64 const > const localSolution,
                                  arrayView1d< real64 const > const pressure,
                                  arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                                  arrayView1d< real64 > pressureScalingFactor,
                                  arrayView1d< real64 > compDensScalingFactor )
    : Base( rankOffset,
            numComp,
            dofKey,
            subRegion,
            localSolution,
            pressure,
            compDens,
            pressureScalingFactor,
            compDensScalingFactor ),
    m_maxRelativePresChange( maxRelativePresChange ),
    m_maxAbsolutePresChange( maxAbsolutePresChange ),
    m_maxCompFracChange( maxCompFracChange ),
    m_maxRelativeCompDensChange( maxRelativeCompDensChange )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables located on the stack
   */
  struct StackVariables : public Base::StackVariables
  {
    GEOS_HOST_DEVICE
    StackVariables()
    { }

    StackVariables( real64 _localMinVal,
                    real64 _localMaxDeltaPres,
                    real64 _localMaxDeltaTemp,
                    real64 _localMaxDeltaCompDens,
                    real64 _localMinPresScalingFactor,
                    real64 _localMinTempScalingFactor,
                    real64 _localMinCompDensScalingFactor )
      :
      Base::StackVariables( _localMinVal ),
      localMaxDeltaPres( _localMaxDeltaPres ),
      localMaxDeltaTemp( _localMaxDeltaTemp ),
      localMaxDeltaCompDens( _localMaxDeltaCompDens ),
      localMinPresScalingFactor( _localMinPresScalingFactor ),
      localMinTempScalingFactor( _localMinTempScalingFactor ),
      localMinCompDensScalingFactor( _localMinCompDensScalingFactor )
    { }

    real64 localMaxDeltaPres;
    real64 localMaxDeltaTemp;
    real64 localMaxDeltaCompDens;

    real64 localMinPresScalingFactor;
    real64 localMinTempScalingFactor;
    real64 localMinCompDensScalingFactor;

  };

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static StackVariables
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > globalScalingFactor( 1.0 );

    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaPres( 0.0 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaTemp( 0.0 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaCompDens( 0.0 );

    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minPresScalingFactor( 1.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minTempScalingFactor( 1.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minCompDensScalingFactor( 1.0 );

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.ghostRank( ei ) >= 0 )
      {
        return;
      }

      StackVariables stack;
      kernelComponent.setup( ei, stack );
      kernelComponent.compute( ei, stack );

      globalScalingFactor.min( stack.localMinVal );

      maxDeltaPres.max( stack.localMaxDeltaPres );
      maxDeltaTemp.max( stack.localMaxDeltaTemp );
      maxDeltaCompDens.max( stack.localMaxDeltaCompDens );

      minPresScalingFactor.min( stack.localMinPresScalingFactor );
      minTempScalingFactor.min( stack.localMinTempScalingFactor );
      minCompDensScalingFactor.min( stack.localMinCompDensScalingFactor );
    } );

    return StackVariables( globalScalingFactor.get(),
                           maxDeltaPres.get(),
                           maxDeltaTemp.get(),
                           maxDeltaCompDens.get(),
                           minPresScalingFactor.get(),
                           minTempScalingFactor.get(),
                           minCompDensScalingFactor.get() );
  }

  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    Base::setup( ei, stack );

    stack.localMaxDeltaPres = 0.0;
    stack.localMaxDeltaTemp = 0.0;
    stack.localMaxDeltaCompDens = 0.0;

    stack.localMinPresScalingFactor = 1.0;
    stack.localMinTempScalingFactor = 1.0;
    stack.localMinCompDensScalingFactor = 1.0;
  }

  /**
   * @brief Compute the local value
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                StackVariables & stack ) const
  {
    computeScalingFactor( ei, stack );
  }

  /**
   * @brief Compute the local value of the scaling factor
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeScalingFactor( localIndex const ei,
                             StackVariables & stack,
                             FUNC && kernelOp = NoOpFunc{} ) const
  {
    real64 constexpr eps = minDensForDivision;

    // compute the change in pressure
    real64 const pres = m_pressure[ei];
    real64 const absPresChange = LvArray::math::abs( m_localSolution[stack.localRow] );
    if( stack.localMaxDeltaPres < absPresChange )
    {
      stack.localMaxDeltaPres = absPresChange;
    }

    // compute pressure scaling factor
    real64 presScalingFactor = 1.0;
    // when enabled, absolute change scaling has a priority over relative change
    if( m_maxAbsolutePresChange > 0.0 ) // maxAbsolutePresChange <= 0.0 means that absolute scaling is disabled
    {
      if( absPresChange > m_maxAbsolutePresChange )
      {
        presScalingFactor = m_maxAbsolutePresChange / absPresChange;
      }
    }
    else if( pres > eps )
    {
      real64 const relativePresChange = absPresChange / pres;
      if( relativePresChange > m_maxRelativePresChange )
      {
        presScalingFactor = m_maxRelativePresChange / relativePresChange;
      }
    }
    m_pressureScalingFactor[ei] = presScalingFactor;
    if( stack.localMinVal > presScalingFactor )
    {
      stack.localMinVal = presScalingFactor;
    }
    if( stack.localMinPresScalingFactor > presScalingFactor )
    {
      stack.localMinPresScalingFactor = presScalingFactor;
    }

    real64 prevTotalDens = 0;
    for( integer ic = 0; ic < m_numComp; ++ic )
    {
      prevTotalDens += m_compDens[ei][ic];
    }

    m_compDensScalingFactor[ei] = 1.0;

    // compute the change in component densities and component fractions
    for( integer ic = 0; ic < m_numComp; ++ic )
    {
      // compute scaling factor based on relative change in component densities
      real64 const absCompDensChange = LvArray::math::abs( m_localSolution[stack.localRow + ic + 1] );
      if( stack.localMaxDeltaCompDens < absCompDensChange )
      {
        stack.localMaxDeltaCompDens = absCompDensChange;
      }

      // This actually checks the change in component fraction, using a lagged total density
      // Indeed we can rewrite the following check as:
      //    | prevCompDens / prevTotalDens - newCompDens / prevTotalDens | > maxCompFracChange
      // Note that the total density in the second term is lagged (i.e, we use prevTotalDens)
      // because I found it more robust than using directly newTotalDens (which can vary also
      // wildly when the compDens change is large)
      real64 const maxAbsCompDensChange = m_maxCompFracChange * prevTotalDens;
      if( absCompDensChange > maxAbsCompDensChange && absCompDensChange > eps )
      {
        real64 const compScalingFactor = maxAbsCompDensChange / absCompDensChange;
        m_compDensScalingFactor[ei] = LvArray::math::min( m_compDensScalingFactor[ei], compScalingFactor );
        if( stack.localMinVal > compScalingFactor )
        {
          stack.localMinVal = compScalingFactor;
        }
        if( stack.localMinCompDensScalingFactor > compScalingFactor )
        {
          stack.localMinCompDensScalingFactor = compScalingFactor;
        }
      }

      // switch from relative to absolute when value is < 1.0
      real64 const maxRelCompDensChange = m_maxRelativeCompDensChange * LvArray::math::max( m_compDens[ei][ic], 1.0 );
      if( absCompDensChange > maxRelCompDensChange && absCompDensChange > eps )
      {
        real64 const compScalingFactor = maxRelCompDensChange / absCompDensChange;
        m_compDensScalingFactor[ei] = LvArray::math::min( m_compDensScalingFactor[ei], compScalingFactor );
        if( stack.localMinVal > compScalingFactor )
        {
          stack.localMinVal = compScalingFactor;
        }
        if( stack.localMinCompDensScalingFactor > compScalingFactor )
        {
          stack.localMinCompDensScalingFactor = compScalingFactor;
        }
      }
    }

    // compute the scaling factor for other vars, such as temperature
    kernelOp();
  }

protected:

  /// Max allowed changes in primary variables
  real64 const m_maxRelativePresChange;
  real64 const m_maxAbsolutePresChange;
  real64 const m_maxCompFracChange;
  real64 const m_maxRelativeCompDensChange;

};

/**
 * @class ScalingForSystemSolutionKernelFactory
 */
class ScalingForSystemSolutionKernelFactory
{
public:

  /*
   * @brief Create and launch the kernel computing the scaling factor
   * @tparam POLICY the kernel policy
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxAbsolutePresChange the max allowed absolute pressure change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] maxRelativeCompDensChange the max allowed comp density change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @return the scaling factor
   */
  template< typename POLICY >
  static ScalingForSystemSolutionKernel::StackVariables
  createAndLaunch( real64 const maxRelativePresChange,
                   real64 const maxAbsolutePresChange,
                   real64 const maxCompFracChange,
                   real64 const maxRelativeCompDensChange,
                   arrayView1d< real64 const > const pressure,
                   arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                   arrayView1d< real64 > pressureScalingFactor,
                   arrayView1d< real64 > compDensScalingFactor,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase & subRegion,
                   arrayView1d< real64 const > const localSolution )
  {

    ScalingForSystemSolutionKernel kernel( maxRelativePresChange, maxAbsolutePresChange, maxCompFracChange, maxRelativeCompDensChange, rankOffset,
                                           numComp, dofKey, subRegion, localSolution, pressure, compDens, pressureScalingFactor, compDensScalingFactor );
    return ScalingForSystemSolutionKernel::launch< POLICY >( subRegion.size(), kernel );
  }
};

/******************************** SolutionCheckKernel ********************************/

/**
 * @class SolutionCheckKernel
 * @brief Define the kernel for checking the updated solution
 */
class SolutionCheckKernel : public ScalingAndCheckingSystemSolutionKernelBase< integer >
{
public:

  using Base = ScalingAndCheckingSystemSolutionKernelBase< integer >;
  using Base::m_rankOffset;
  using Base::m_numComp;
  using Base::m_dofNumber;
  using Base::m_ghostRank;
  using Base::m_localSolution;
  using Base::m_pressure;
  using Base::m_compDens;

  /**
   * @brief Create a new kernel instance
   * @param[in] allowCompDensChopping flag to allow the component density chopping
   * @param[in] scalingFactor the scaling factor
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @param[in] pressure the pressure vector
   * @param[in] compDens the component density vector
   */
  SolutionCheckKernel( integer const allowCompDensChopping,
                       integer const allowNegativePressure,
                       CompositionalMultiphaseFVM::ScalingType const scalingType,
                       real64 const scalingFactor,
                       arrayView1d< real64 const > const pressure,
                       arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                       arrayView1d< real64 > pressureScalingFactor,
                       arrayView1d< real64 > compDensScalingFactor,
                       globalIndex const rankOffset,
                       integer const numComp,
                       string const dofKey,
                       ElementSubRegionBase const & subRegion,
                       arrayView1d< real64 const > const localSolution )
    : Base( rankOffset,
            numComp,
            dofKey,
            subRegion,
            localSolution,
            pressure,
            compDens,
            pressureScalingFactor,
            compDensScalingFactor ),
    m_allowCompDensChopping( allowCompDensChopping ),
    m_allowNegativePressure( allowNegativePressure ),
    m_scalingFactor( scalingFactor ),
    m_scalingType( scalingType )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables located on the stack
   */
  struct StackVariables : public Base::StackVariables
  {
    GEOS_HOST_DEVICE
    StackVariables()
    { }

    StackVariables( real64 _localMinVal,
                    real64 _localMinPres,
                    real64 _localMinDens,
                    real64 _localMinTotalDens,
                    integer _localNumNegPressures,
                    integer _localNumNegDens,
                    integer _localNumNegTotalDens )
      :
      Base::StackVariables( _localMinVal ),
      localMinPres( _localMinPres ),
      localMinDens( _localMinDens ),
      localMinTotalDens( _localMinTotalDens ),
      localNumNegPressures( _localNumNegPressures ),
      localNumNegDens( _localNumNegDens ),
      localNumNegTotalDens( _localNumNegTotalDens )
    { }

    real64 localMinPres;
    real64 localMinDens;
    real64 localMinTotalDens;

    integer localNumNegPressures;
    integer localNumNegDens;
    integer localNumNegTotalDens;

  };

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static StackVariables
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, integer > globalMinVal( 1 );

    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minPres( 0.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minDens( 0.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minTotalDens( 0.0 );

    RAJA::ReduceSum< ReducePolicy< POLICY >, integer > numNegPressures( 0 );
    RAJA::ReduceSum< ReducePolicy< POLICY >, integer > numNegDens( 0 );
    RAJA::ReduceSum< ReducePolicy< POLICY >, integer > numNegTotalDens( 0 );

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.ghostRank( ei ) >= 0 )
      {
        return;
      }

      StackVariables stack;
      kernelComponent.setup( ei, stack );
      kernelComponent.compute( ei, stack );

      globalMinVal.min( stack.localMinVal );

      minPres.min( stack.localMinPres );
      minDens.min( stack.localMinDens );
      minTotalDens.min( stack.localMinTotalDens );

      numNegPressures += stack.localNumNegPressures;
      numNegDens += stack.localNumNegDens;
      numNegTotalDens += stack.localNumNegTotalDens;
    } );

    return StackVariables( globalMinVal.get(),
                           minPres.get(),
                           minDens.get(),
                           minTotalDens.get(),
                           numNegPressures.get(),
                           numNegDens.get(),
                           numNegTotalDens.get() );
  }

  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    Base::setup( ei, stack );

    stack.localMinPres = 0.0;
    stack.localMinDens = 0.0;
    stack.localMinTotalDens = 0.0;

    stack.localNumNegPressures = 0;
    stack.localNumNegDens = 0;
    stack.localNumNegTotalDens = 0;
  }

  /**
   * @brief Compute the local value
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                StackVariables & stack ) const
  {
    computeSolutionCheck( ei, stack );
  }

  /**
   * @brief Compute the local value of the check
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeSolutionCheck( localIndex const ei,
                             StackVariables & stack,
                             FUNC && kernelOp = NoOpFunc{} ) const
  {
    bool const localScaling = m_scalingType == CompositionalMultiphaseFVM::ScalingType::Local;

    real64 const newPres = m_pressure[ei] + (localScaling ? m_pressureScalingFactor[ei] : m_scalingFactor) * m_localSolution[stack.localRow];
    if( newPres < 0 )
    {
      if( !m_allowNegativePressure )
      {
        stack.localMinVal = 0;
      }
      stack.localNumNegPressures += 1;
      if( newPres < stack.localMinPres )
        stack.localMinPres = newPres;
    }

    // if component density chopping is not allowed, the time step fails if a component density is negative
    // otherwise, we just check that the total density is positive, and negative component densities
    // will be chopped (i.e., set to zero) in ApplySystemSolution)
    if( !m_allowCompDensChopping )
    {
      for( integer ic = 0; ic < m_numComp; ++ic )
      {
        real64 const newDens = m_compDens[ei][ic] + (localScaling ? m_compDensScalingFactor[ei] : m_scalingFactor) * m_localSolution[stack.localRow + ic + 1];
        if( newDens < 0 )
        {
          stack.localMinVal = 0;
          stack.localNumNegDens += 1;
          if( newDens < stack.localMinDens )
            stack.localMinDens = newDens;
        }
      }
    }
    else
    {
      real64 totalDens = 0.0;
      for( integer ic = 0; ic < m_numComp; ++ic )
      {
        real64 const newDens = m_compDens[ei][ic] + (localScaling ? m_compDensScalingFactor[ei] : m_scalingFactor) * m_localSolution[stack.localRow + ic + 1];
        totalDens += ( newDens > 0.0 ) ? newDens : 0.0;
      }
      if( totalDens < 0 )
      {
        stack.localMinVal = 0;
        stack.localNumNegTotalDens += 1;
        if( totalDens < stack.localMinTotalDens )
          stack.localMinTotalDens = totalDens;
      }
    }

    kernelOp();
  }

protected:

  /// flag to allow the component density chopping
  integer const m_allowCompDensChopping;

  /// flag to allow negative pressure values
  integer const m_allowNegativePressure;

  /// scaling factor
  real64 const m_scalingFactor;

  /// scaling type (global or local)
  CompositionalMultiphaseFVM::ScalingType const m_scalingType;

};

/**
 * @class SolutionCheckKernelFactory
 */
class SolutionCheckKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] allowCompDensChopping flag to allow the component density chopping
   * @param[in] scalingFactor the scaling factor
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   */
  template< typename POLICY >
  static SolutionCheckKernel::StackVariables
  createAndLaunch( integer const allowCompDensChopping,
                   integer const allowNegativePressure,
                   CompositionalMultiphaseFVM::ScalingType const scalingType,
                   real64 const scalingFactor,
                   arrayView1d< real64 const > const pressure,
                   arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                   arrayView1d< real64 > pressureScalingFactor,
                   arrayView1d< real64 > compDensScalingFactor,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase & subRegion,
                   arrayView1d< real64 const > const localSolution )
  {

    SolutionCheckKernel kernel( allowCompDensChopping, allowNegativePressure, scalingType, scalingFactor,
                                pressure, compDens, pressureScalingFactor, compDensScalingFactor, rankOffset,
                                numComp, dofKey, subRegion, localSolution );
    return SolutionCheckKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< 2 >
{
public:

  using Base = physicsSolverBaseKernels::ResidualNormKernelBase< 2 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      integer const numComponents,
                      ElementSubRegionBase const & subRegion,
                      constitutive::MultiFluidBase const & fluid,
                      constitutive::CoupledSolidBase const & solid,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_numComponents( numComponents ),
    m_volume( subRegion.getElementVolume() ),
    m_porosity_n( solid.getPorosity_n() ),
    m_totalDens_n( fluid.totalDensity_n() )
  {}

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    // this should never be zero if the simulation is set up correctly, but we never know
    real64 const massNormalizer = LvArray::math::max( m_minNormalizer, m_totalDens_n[ei][0] * m_porosity_n[ei][0] * m_volume[ei] );
    real64 const volumeNormalizer = LvArray::math::max( m_minNormalizer, m_porosity_n[ei][0] * m_volume[ei] );

    // step 1: mass residuals

    for( integer idof = 0; idof < m_numComponents; ++idof )
    {
      real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / massNormalizer;
      if( valMass > stack.localValue[0] )
      {
        stack.localValue[0] = valMass;
      }
    }

    // step 2: volume residual

    real64 const valVol = LvArray::math::abs( m_localResidual[stack.localRow + m_numComponents] ) / volumeNormalizer;
    if( valVol > stack.localValue[1] )
    {
      stack.localValue[1] = valVol;
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    // note: for the L2 norm, we bundle the volume and mass residuals/normalizers

    real64 const massNormalizer = LvArray::math::max( m_minNormalizer, m_totalDens_n[ei][0] * m_porosity_n[ei][0] * m_volume[ei] );

    // step 1: mass residuals

    for( integer idof = 0; idof < m_numComponents; ++idof )
    {
      stack.localValue[0] += m_localResidual[stack.localRow + idof] * m_localResidual[stack.localRow + idof];
      stack.localNormalizer[0] += massNormalizer;
    }

    // step 2: volume residual

    real64 const val = m_localResidual[stack.localRow + m_numComponents] * m_totalDens_n[ei][0]; // we need a mass here, hence the
                                                                                                 // multiplication
    stack.localValue[1] += val * val;
    stack.localNormalizer[1] += massNormalizer;
  }


protected:

  /// Number of fluid coponents
  integer const m_numComponents;

  /// View on the volume
  arrayView1d< real64 const > const m_volume;

  /// View on porosity at the previous converged time step
  arrayView2d< real64 const > const m_porosity_n;

  /// View on total mass/molar density at the previous converged time step
  arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const m_totalDens_n;

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
   * @param[in] numComps the number of fluid components
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
                   integer const numComps,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   constitutive::MultiFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[2],
                   real64 (& residualNormalizer)[2] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank, numComps, subRegion, fluid, solid, minNormalizer );
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

/******************************** StatisticsKernel ********************************/

struct StatisticsKernel
{
  template< typename POLICY >
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

  template< typename POLICY >
  static void
  launch( localIndex const size,
          integer const numComps,
          integer const numPhases,
          real64 const relpermThreshold,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & deltaPres,
          arrayView1d< real64 const > const & temp,
          arrayView1d< real64 const > const & refPorosity,
          arrayView2d< real64 const > const & porosity,
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDensity,
          arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > const & phaseCompFraction,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > const & phaseTrappedVolFrac,
          arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > const & phaseRelperm,
          real64 & minPres,
          real64 & avgPresNumerator,
          real64 & maxPres,
          real64 & minDeltaPres,
          real64 & maxDeltaPres,
          real64 & minTemp,
          real64 & avgTempNumerator,
          real64 & maxTemp,
          real64 & totalUncompactedPoreVol,
          arrayView1d< real64 > const & phaseDynamicPoreVol,
          arrayView1d< real64 > const & phaseMass,
          arrayView1d< real64 > const & trappedPhaseMass,
          arrayView1d< real64 > const & immobilePhaseMass,
          arrayView2d< real64 > const & dissolvedComponentMass )
  {
    RAJA::ReduceMin< parallelDeviceReduce, real64 > subRegionMinPres( LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionAvgPresNumerator( 0.0 );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxPres( -LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real64 > subRegionMinDeltaPres( LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxDeltaPres( -LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real64 > subRegionMinTemp( LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionAvgTempNumerator( 0.0 );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxTemp( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionTotalUncompactedPoreVol( 0.0 );

    // For this arrays phaseDynamicPoreVol, phaseMass, dissolvedComponentMass,
    // using an array of ReduceSum leads to a formal parameter overflow in CUDA.
    // As a workaround, we use a slice with RAJA::atomicAdd instead

    forAll< parallelDevicePolicy<> >( size, [numComps,
                                             numPhases,
                                             relpermThreshold,
                                             elemGhostRank,
                                             volume,
                                             refPorosity,
                                             porosity,
                                             pres,
                                             deltaPres,
                                             temp,
                                             phaseDensity,
                                             phaseVolFrac,
                                             phaseTrappedVolFrac,
                                             phaseRelperm,
                                             phaseCompFraction,
                                             subRegionMinPres,
                                             subRegionAvgPresNumerator,
                                             subRegionMaxPres,
                                             subRegionMinDeltaPres,
                                             subRegionMaxDeltaPres,
                                             subRegionMinTemp,
                                             subRegionAvgTempNumerator,
                                             subRegionMaxTemp,
                                             subRegionTotalUncompactedPoreVol,
                                             phaseDynamicPoreVol,
                                             phaseMass,
                                             trappedPhaseMass,
                                             immobilePhaseMass,
                                             dissolvedComponentMass] GEOS_HOST_DEVICE ( localIndex const ei )
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

      subRegionMaxDeltaPres.max( deltaPres[ei] );
      subRegionMinDeltaPres.min( deltaPres[ei] );

      subRegionMinTemp.min( temp[ei] );
      subRegionAvgTempNumerator += uncompactedPoreVol * temp[ei];
      subRegionMaxTemp.max( temp[ei] );
      subRegionTotalUncompactedPoreVol += uncompactedPoreVol;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        real64 const elemPhaseVolume = dynamicPoreVol * phaseVolFrac[ei][ip];
        real64 const elemPhaseMass = phaseDensity[ei][0][ip] * elemPhaseVolume;
        real64 const elemTrappedPhaseMass = phaseDensity[ei][0][ip] * dynamicPoreVol * phaseTrappedVolFrac[ei][0][ip];
        // RAJA::atomicAdd used here because we do not use ReduceSum here (for the reason explained above)
        RAJA::atomicAdd( parallelDeviceAtomic{}, &phaseDynamicPoreVol[ip], elemPhaseVolume );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &phaseMass[ip], elemPhaseMass );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &trappedPhaseMass[ip], elemTrappedPhaseMass );
        if( phaseRelperm[ei][0][ip] < relpermThreshold )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &immobilePhaseMass[ip], elemPhaseMass );
        }
        for( integer ic = 0; ic < numComps; ++ic )
        {
          // RAJA::atomicAdd used here because we do not use ReduceSum here (for the reason explained above)
          RAJA::atomicAdd( parallelDeviceAtomic{}, &dissolvedComponentMass[ip][ic], phaseCompFraction[ei][0][ip][ic] * elemPhaseMass );
        }
      }

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

    // dummy loop to bring data back to the CPU
    forAll< serialPolicy >( 1, [phaseDynamicPoreVol, phaseMass, trappedPhaseMass, immobilePhaseMass, dissolvedComponentMass] ( localIndex const )
    {
      GEOS_UNUSED_VAR( phaseDynamicPoreVol, phaseMass, trappedPhaseMass, immobilePhaseMass, dissolvedComponentMass );
    } );
  }
};

/******************************** HydrostaticPressureKernel ********************************/

struct HydrostaticPressureKernel
{

  // TODO: this type of constants should be centralized somewhere or provided by fluid model
  static real64 constexpr MIN_FOR_PHASE_PRESENCE = 1e-12;

  enum class ReturnType : integer
  {
    FAILED_TO_CONVERGE = 0,
    DETECTED_MULTIPHASE_FLOW = 1,
    SUCCESS = 2
  };

  template< typename FLUID_WRAPPER >
  static ReturnType
  computeHydrostaticPressure( integer const numComps,
                              integer const numPhases,
                              integer const ipInit,
                              integer const maxNumEquilIterations,
                              real64 const & equilTolerance,
                              real64 const (&gravVector)[ 3 ],
                              FLUID_WRAPPER fluidWrapper,
                              arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                              TableFunction::KernelWrapper tempTableWrapper,
                              real64 const & refElevation,
                              real64 const & refPres,
                              arraySlice1d< real64 const > const & refPhaseMassDens,
                              real64 const & newElevation,
                              real64 & newPres,
                              arraySlice1d< real64 > const & newPhaseMassDens )
  {
    // fluid properties at this elevation
    StackArray< real64, 2, constitutive::MultiFluidBase::MAX_NUM_COMPONENTS, compflow::LAYOUT_COMP > compFrac( 1, numComps );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > phaseFrac( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > phaseDens( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > phaseMassDens( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > phaseVisc( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > phaseEnthalpy( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > phaseInternalEnergy( 1, 1, numPhases );
    StackArray< real64, 4, constitutive::MultiFluidBase::MAX_NUM_PHASES *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS,
                constitutive::multifluid::LAYOUT_PHASE_COMP > phaseCompFrac( 1, 1, numPhases, numComps );
    real64 totalDens = 0.0;

    bool isSinglePhaseFlow = true;

    // Step 1: compute the hydrostatic pressure at the current elevation

    real64 const gravCoef = gravVector[2] * ( refElevation - newElevation );
    real64 const temp = tempTableWrapper.compute( &newElevation );
    for( integer ic = 0; ic < numComps; ++ic )
    {
      compFrac[0][ic] = compFracTableWrappers[ic].compute( &newElevation );
    }

    // Step 2: guess the pressure with the refPhaseMassDensity

    real64 pres0 = refPres - refPhaseMassDens[ipInit] * gravCoef;
    real64 pres1 = 0.0;

    // Step 3: compute the mass density at this elevation using the guess, and update pressure

    constitutive::MultiFluidBase::KernelWrapper::computeValues( fluidWrapper,
                                                                pres0,
                                                                temp,
                                                                compFrac[0],
                                                                phaseFrac[0][0],
                                                                phaseDens[0][0],
                                                                phaseMassDens[0][0],
                                                                phaseVisc[0][0],
                                                                phaseEnthalpy[0][0],
                                                                phaseInternalEnergy[0][0],
                                                                phaseCompFrac[0][0],
                                                                totalDens );
    pres1 = refPres - 0.5 * ( refPhaseMassDens[ipInit] + phaseMassDens[0][0][ipInit] ) * gravCoef;

    // Step 4: fixed-point iteration until convergence

    bool equilHasConverged = false;
    for( integer eqIter = 0; eqIter < maxNumEquilIterations; ++eqIter )
    {

      // check convergence
      equilHasConverged = ( LvArray::math::abs( pres0 - pres1 ) < equilTolerance );
      pres0 = pres1;

      // if converged, check number of phases and move on
      if( equilHasConverged )
      {
        // make sure that the fluid is single-phase, other we have to issue a warning (for now)
        // if only one phase is mobile, we are in good shape (unfortunately it is hard to access relperm from here)
        localIndex numberOfPhases = 0;
        for( integer ip = 0; ip < numPhases; ++ip )
        {
          if( phaseFrac[0][0][ip] > MIN_FOR_PHASE_PRESENCE )
          {
            numberOfPhases++;
          }
        }
        if( numberOfPhases > 1 )
        {
          isSinglePhaseFlow = false;
        }

        break;
      }

      // compute the mass density at this elevation using the previous pressure, and compute the new pressure
      constitutive::MultiFluidBase::KernelWrapper::computeValues( fluidWrapper,
                                                                  pres0,
                                                                  temp,
                                                                  compFrac[0],
                                                                  phaseFrac[0][0],
                                                                  phaseDens[0][0],
                                                                  phaseMassDens[0][0],
                                                                  phaseVisc[0][0],
                                                                  phaseEnthalpy[0][0],
                                                                  phaseInternalEnergy[0][0],
                                                                  phaseCompFrac[0][0],
                                                                  totalDens );
      pres1 = refPres - 0.5 * ( refPhaseMassDens[ipInit] + phaseMassDens[0][0][ipInit] ) * gravCoef;
    }

    // Step 5: save the hydrostatic pressure and the corresponding density

    newPres = pres1;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      newPhaseMassDens[ip] = phaseMassDens[0][0][ip];
    }

    if( !equilHasConverged )
    {
      return ReturnType::FAILED_TO_CONVERGE;
    }
    else if( !isSinglePhaseFlow )
    {
      return ReturnType::DETECTED_MULTIPHASE_FLOW;
    }
    else
    {
      return ReturnType::SUCCESS;
    }
  }

  template< typename FLUID_WRAPPER >
  static ReturnType
  launch( localIndex const size,
          integer const numComps,
          integer const numPhases,
          integer const ipInit,
          integer const maxNumEquilIterations,
          real64 const equilTolerance,
          real64 const (&gravVector)[ 3 ],
          real64 const & minElevation,
          real64 const & elevationIncrement,
          real64 const & datumElevation,
          real64 const & datumPres,
          FLUID_WRAPPER fluidWrapper,
          arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
          TableFunction::KernelWrapper tempTableWrapper,
          arrayView1d< arrayView1d< real64 > const > elevationValues,
          arrayView1d< real64 > pressureValues )
  {

    ReturnType returnVal = ReturnType::SUCCESS;

    // Step 1: compute the phase mass densities at datum

    // datum fluid properties
    array2d< real64, compflow::LAYOUT_COMP > datumCompFrac( 1, numComps );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPhaseFrac( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPhaseDens( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPhaseMassDens( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPhaseVisc( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPhaseEnthalpy( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPhaseInternalEnergy( 1, 1, numPhases );
    array4d< real64, constitutive::multifluid::LAYOUT_PHASE_COMP > datumPhaseCompFrac( 1, 1, numPhases, numComps );
    real64 datumTotalDens = 0.0;

    real64 const datumTemp = tempTableWrapper.compute( &datumElevation );
    for( integer ic = 0; ic < numComps; ++ic )
    {
      datumCompFrac[0][ic] = compFracTableWrappers[ic].compute( &datumElevation );
    }
    constitutive::MultiFluidBase::KernelWrapper::computeValues( fluidWrapper,
                                                                datumPres,
                                                                datumTemp,
                                                                datumCompFrac[0],
                                                                datumPhaseFrac[0][0],
                                                                datumPhaseDens[0][0],
                                                                datumPhaseMassDens[0][0],
                                                                datumPhaseVisc[0][0],
                                                                datumPhaseEnthalpy[0][0],
                                                                datumPhaseInternalEnergy[0][0],
                                                                datumPhaseCompFrac[0][0],
                                                                datumTotalDens );

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

    array2d< real64 > phaseMassDens( pressureValues.size(), numPhases );
    // temporary array without permutation to compile on Lassen
    array1d< real64 > datumPhaseMassDensTmp( numPhases );
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      datumPhaseMassDensTmp[ip] = datumPhaseMassDens[0][0][ip];
    }

    ReturnType const refReturnVal =
      computeHydrostaticPressure( numComps,
                                  numPhases,
                                  ipInit,
                                  maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluidWrapper,
                                  compFracTableWrappers,
                                  tempTableWrapper,
                                  datumElevation,
                                  datumPres,
                                  datumPhaseMassDensTmp,
                                  elevationValues[0][iRef],
                                  pressureValues[iRef],
                                  phaseMassDens[iRef] );
    if( refReturnVal == ReturnType::FAILED_TO_CONVERGE )
    {
      return ReturnType::FAILED_TO_CONVERGE;
    }
    else if( refReturnVal == ReturnType::DETECTED_MULTIPHASE_FLOW )
    {
      returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
    }

    // Step 4: for each elevation above the reference elevation, compute the pressure

    localIndex const numEntriesAboveRef = size - iRef - 1;
    forAll< serialPolicy >( numEntriesAboveRef, [=, &returnVal] ( localIndex const i )
    {
      ReturnType const returnValAboveRef =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ipInit,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    elevationValues[0][iRef+i],
                                    pressureValues[iRef+i],
                                    phaseMassDens[iRef+i],
                                    elevationValues[0][iRef+i+1],
                                    pressureValues[iRef+i+1],
                                    phaseMassDens[iRef+i+1] );
      if( returnValAboveRef == ReturnType::FAILED_TO_CONVERGE )
      {
        returnVal = ReturnType::FAILED_TO_CONVERGE;
      }
      else if( ( returnValAboveRef == ReturnType::DETECTED_MULTIPHASE_FLOW ) &&
               ( returnVal != ReturnType::FAILED_TO_CONVERGE ) )
      {
        returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
      }

    } );

    // Step 5: for each elevation below the reference elevation, compute the pressure

    localIndex const numEntriesBelowRef = iRef;
    forAll< serialPolicy >( numEntriesBelowRef, [=, &returnVal] ( localIndex const i )
    {
      ReturnType const returnValBelowRef =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ipInit,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    elevationValues[0][iRef-i],
                                    pressureValues[iRef-i],
                                    phaseMassDens[iRef-i],
                                    elevationValues[0][iRef-i-1],
                                    pressureValues[iRef-i-1],
                                    phaseMassDens[iRef-i-1] );
      if( returnValBelowRef == ReturnType::FAILED_TO_CONVERGE )
      {
        returnVal = ReturnType::FAILED_TO_CONVERGE;
      }
      else if( ( returnValBelowRef == ReturnType::DETECTED_MULTIPHASE_FLOW ) &&
               ( returnVal != ReturnType::FAILED_TO_CONVERGE ) )
      {
        returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
      }

    } );

    return returnVal;
  }

};


/******************************** Kernel launch machinery ********************************/


template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelectorCompTherm( integer const numComp, bool const isThermal, ARGS && ... args )
{
  geos::internal::kernelLaunchSelectorCompThermSwitch( numComp, isThermal, [&] ( auto NC, auto ISTHERMAL )
  {
    KERNELWRAPPER::template launch< NC(), ISTHERMAL() >( std::forward< ARGS >( args )... );
  } );
}

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector1( integer const numComp, ARGS && ... args )
{
  internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
  {
    KERNELWRAPPER::template launch< NC() >( std::forward< ARGS >( args )... );
  } );
}

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector2( integer const numComp, integer const numPhase, ARGS && ... args )
{
  // Ideally this would be inside the dispatch, but it breaks on Summit with GCC 9.1.0 and CUDA 11.0.3.
  if( numPhase == 2 )
  {
    internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
    {
      KERNELWRAPPER::template launch< NC(), 2 >( std::forward< ARGS >( args ) ... );
    } );
  }
  else if( numPhase == 3 )
  {
    internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
    {
      KERNELWRAPPER::template launch< NC(), 3 >( std::forward< ARGS >( args ) ... );
    } );
  }
  else
  {
    GEOS_ERROR( "Unsupported number of phases: " << numPhase );
  }
}

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector_NC_NP_THERM( integer const numComp, integer const numPhase, integer const isThermal, ARGS && ... args )
{
  // Ideally this would be inside the dispatch, but it breaks on Summit with GCC 9.1.0 and CUDA 11.0.3.
  if( isThermal )
  {
    if( numPhase == 2 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        KERNELWRAPPER::template launch< NC(), 2, 1 >( std::forward< ARGS >( args ) ... );
      } );
    }
    else if( numPhase == 3 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        KERNELWRAPPER::template launch< NC(), 3, 1 >( std::forward< ARGS >( args ) ... );
      } );
    }
    else
    {
      GEOS_ERROR( "Unsupported number of phases: " << numPhase );
    }
  }
  else
  {
    if( numPhase == 2 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        KERNELWRAPPER::template launch< NC(), 2, 0 >( std::forward< ARGS >( args ) ... );
      } );
    }
    else if( numPhase == 3 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        KERNELWRAPPER::template launch< NC(), 3, 0 >( std::forward< ARGS >( args ) ... );
      } );
    }
    else
    {
      GEOS_ERROR( "Unsupported number of phases: " << numPhase );
    }
  }
}

} // namespace isothermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEBASEKERNELS_HPP
