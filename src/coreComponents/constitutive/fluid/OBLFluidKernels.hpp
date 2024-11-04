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
 * @file OBLFluidKernels.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_OBLFLUIDKERNELS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_OBLFLUIDKERNELS_HPP_

#include "functions/MultilinearInterpolatorStaticKernels.hpp"
#include "functions/MultilinearInterpolatorAdaptiveKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/ReactiveCompositionalMultiphaseOBLFields.hpp"

namespace geos
{

namespace oblFluidKernels
{

using namespace constitutive;

// The number of operators in use depends on:
// 1. Number of phases
// 2. Number of components
// 3. Features required in simulation (now its only one - Energy balance)
// This number needs to be used in solver and in kernels (as a template parameter)
// IMHO, this number is too big ( order of 10-100) to be treated via a kernelLaunchSelectorSwitch construct
// Hence, a way to define it once and for all is needed.
// Could be constexpr member of solver, but passing constexpr lambdas require -std=c++17 if I`m not mistaken

constexpr integer COMPUTE_NUM_OPS ( integer const NP, integer const NC, bool ENERGY )
{
  auto DOF = NC + ENERGY;
  return DOF /*accumulation*/ +
         DOF * NP /*flux*/ +
         NP /*up_constant*/ +
         DOF * NP /*gradient*/ +
         DOF /*kinetic rate*/ +
         2 /*rock internal energy and conduction*/ +
         2 * NP /*gravity and capillarity*/ +
         1 /*rock porosity*/ +
         1;
}


/******************************** Kernel launch machinery ********************************/
namespace internal
{

template< bool ENABLE_ENERGY, integer NUM_PHASES, typename T, typename LAMBDA >
void kernelLaunchSelectorCompSwitch( T numComps, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorCompSwitch: type should be integral" );

  switch( numComps )
  {
    case 1:
    { lambda( std::integral_constant< T, NUM_PHASES >(), std::integral_constant< T, 1 >(), std::integral_constant< bool, ENABLE_ENERGY >() ); return; }
    case 2:
    { lambda( std::integral_constant< T, NUM_PHASES >(), std::integral_constant< T, 2 >(), std::integral_constant< bool, ENABLE_ENERGY >() ); return; }
    case 3:
    { lambda( std::integral_constant< T, NUM_PHASES >(), std::integral_constant< T, 3 >(), std::integral_constant< bool, ENABLE_ENERGY >() ); return; }
    case 4:
    { lambda( std::integral_constant< T, NUM_PHASES >(), std::integral_constant< T, 4 >(), std::integral_constant< bool, ENABLE_ENERGY >() ); return; }
    case 5:
    { lambda( std::integral_constant< T, NUM_PHASES >(), std::integral_constant< T, 5 >(), std::integral_constant< bool, ENABLE_ENERGY >()); return; }
    default:
    { GEOS_ERROR( "Unsupported number of components: " << numComps ); }
  }
}

template< bool ENABLE_ENERGY, typename T, typename LAMBDA >
void kernelLaunchSelectorPhaseSwitch( T numPhases, T numComps, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorPhaseSwitch: type should be integral" );

  switch( numPhases )
  {
    case 1:
    { kernelLaunchSelectorCompSwitch< ENABLE_ENERGY, 1 >( numComps, lambda ); return; }
    case 2:
    { kernelLaunchSelectorCompSwitch< ENABLE_ENERGY, 2 >( numComps, lambda ); return; }
    case 3:
    { kernelLaunchSelectorCompSwitch< ENABLE_ENERGY, 3 >( numComps, lambda ); return; }
    default:
    { GEOS_ERROR( "Unsupported number of phases: " << numPhases ); }
  }
}

template< typename T, typename LAMBDA >
void kernelLaunchSelectorEnergySwitch( T numPhases, T numComps, bool enableEnergyBalance, LAMBDA && lambda )
{
  if( enableEnergyBalance )
  {
    kernelLaunchSelectorPhaseSwitch< true >( numPhases, numComps, lambda );
  }
  else
  {
    kernelLaunchSelectorPhaseSwitch< false >( numPhases, numComps, lambda );
  }
}

} // namespace internal

/******************************** OBLOperatorsKernel ********************************/

/**
 * @class OBLOperatorsKernel
 * @tparam NUM_PHASES number of phases
 * @tparam NUM_COMPS number of components
 * @tparam ENABLE_ENERGY flag if energy balance equation is assembled
 * @brief Compute OBL Operators and derivatives
 */
template< integer NUM_PHASES, integer NUM_COMPS, bool ENABLE_ENERGY >
class OBLOperatorsKernel
{
public:
  /// Compile time value for the energy balance switch
  static constexpr integer enableEnergyBalance = ENABLE_ENERGY;
  /// Compile time value for the number of components
  static constexpr integer numComps = NUM_COMPS;

  /// Compile time value for the number of dimensions
  static constexpr integer numDofs = numComps + enableEnergyBalance;
  // /// Compile time value for the number of operators
  static constexpr integer numOps = COMPUTE_NUM_OPS( NUM_PHASES, NUM_COMPS, ENABLE_ENERGY );

  static constexpr real64 barToPascalMult = 1e5;
  static constexpr real64 pascalToBarMult = 1.0 / 1e5;

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
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] OBLOperatorsTable the OBL table function kernel
   */
  OBLOperatorsKernel( ObjectManagerBase & subRegion,
                      MultilinearInterpolatorBaseKernel< numDofs, numOps > const * OBLOperatorsTable )
    :
    m_OBLOperatorsTable( OBLOperatorsTable ),
    m_pressure( subRegion.getField< fields::flow::pressure >() ),
    m_compFrac( subRegion.getField< fields::flow::globalCompFraction >() ),
    m_temperature( subRegion.getField< fields::flow::temperature >() ),
    m_OBLOperatorValues ( subRegion.getField< fields::flow::OBLOperatorValues >()),
    m_OBLOperatorDerivatives ( subRegion.getField< fields::flow::OBLOperatorDerivatives >())
  {}

  /**
   * @brief Compute the operator values and derivatives for an element
   * @param[in] ei the element index
   */
  GEOS_HOST_DEVICE
  inline
  void compute( localIndex const ei ) const
  {
    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const compFrac = m_compFrac[ei];
    arraySlice1d< real64, compflow::USD_OBL_VAL - 1 > const & OBLVals = m_OBLOperatorValues[ei];
    arraySlice2d< real64, compflow::USD_OBL_DER - 1 > const & OBLDers = m_OBLOperatorDerivatives[ei];
    real64 state[numDofs];

    // we need to convert pressure from Pa (internal unit in GEOSX) to bar (internal unit in DARTS)
    state[0] = m_pressure[ei] * pascalToBarMult;

    // the last component fraction is not used to define the state
    for( integer i = 1; i < numComps; ++i )
    {
      state[i] = compFrac[i - 1];
    }

    if( enableEnergyBalance )
    {
      state[numDofs - 1] = m_temperature[ei];
    }

    m_OBLOperatorsTable->compute( state, OBLVals, OBLDers );

    // we do not perform derivatives unit conversion here:
    // instead we postpone it till all the derivatives are fully formed, and only then apply the factor only once in 'complete' function
    // scaling the whole system might be even better solution (every pressure column needs to be multiplied by pascalToBarMult)
  }

private:

  // inputs
  MultilinearInterpolatorBaseKernel< numDofs, numOps > const * m_OBLOperatorsTable;

  // Views on primary variables and their updates
  arrayView1d< real64 const > m_pressure;
  arrayView2d< real64 const, compflow::USD_COMP > m_compFrac;
  arrayView1d< real64 const > m_temperature;

  // outputs

  // Views on OBL operator values and derivatives
  arrayView2d< real64, compflow::USD_OBL_VAL > m_OBLOperatorValues;
  arrayView3d< real64, compflow::USD_OBL_DER > m_OBLOperatorDerivatives;
};

/**
 * @class OBLOperatorsKernelFactory
 */
class OBLOperatorsKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numPhases the number of phases
   * @param[in] numComponents the number of components
   * @param[in] enableEnergyBalance flag if energy balance equation is assembled
   * @param[in] subRegion the element subregion
   * @param[in] function the OBL table function
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numPhases,
                   integer const numComponents,
                   bool const enableEnergyBalance,
                   ObjectManagerBase & subRegion,
                   constitutive::OBLFluid * oblFluid )
  {
    internal::kernelLaunchSelectorEnergySwitch( numPhases, numComponents, enableEnergyBalance, [&] ( auto NP, auto NC, auto E )
    {
      integer constexpr ENABLE_ENERGY = E();
      integer constexpr NUM_PHASES = NP();
      integer constexpr NUM_COMPS = NC();
      integer constexpr NUM_DIMS = ENABLE_ENERGY + NUM_COMPS;
      integer constexpr NUM_OPS  = COMPUTE_NUM_OPS( NUM_PHASES, NUM_COMPS, ENABLE_ENERGY );

      if( oblFluid->getInterpolatorMode() == constitutive::OBLInterpolatorMode::Static )
      {
        MultivariableTableFunction const & function = oblFluid->getTable();

        GEOS_THROW_IF_NE_MSG( NUM_DIMS, function.numDims(),
                              GEOS_FMT( "The number of degrees of freedom per cell used in the solver has a value of {}, "
                                        "whereas it as a value of {} in the operator table (at {}).",
                                        NUM_DIMS, function.numDims(),
                                        function.getName() ),
                              InputError );

        GEOS_THROW_IF_NE_MSG( NUM_OPS, function.numOps(),
                              GEOS_FMT( "The number of operators per cell used in the solver has a value of {}, "
                                        "whereas it as a value of {} in the operator table (at {}).",
                                        NUM_OPS, function.numOps(),
                                        function.getName() ),
                              InputError );

        MultilinearInterpolatorStaticKernel< NUM_DIMS, NUM_OPS > const interpolationKernel(
          function.getAxisMinimums(),
          function.getAxisMaximums(),
          function.getAxisPoints(),
          function.getAxisSteps(),
          function.getAxisStepInvs(),
          function.getAxisHypercubeMults(),
          function.getHypercubeData()
          );

        OBLOperatorsKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY > kernel( subRegion, &interpolationKernel );
        OBLOperatorsKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY >::template launch< POLICY >( subRegion.size(), kernel );
      }
      else /* if ( oblFluid->getInterpolatorMode() == constitutive::OBLInterpolatorMode::Adaptive ) */
      {
        // Check if Python function was assigned to wrapper
        auto pyFunctionPtr = oblFluid->getPythonFunction();

        GEOS_ERROR_IF( pyFunctionPtr->py_evaluate_func == nullptr,
                       GEOS_FMT( "{}: py_evaluate_func is not specified",
                                 pyFunctionPtr->getName())
                       );
        MultilinearInterpolatorAdaptiveKernel< NUM_DIMS, NUM_OPS > const interpolationKernel( pyFunctionPtr );

        OBLOperatorsKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY > kernel( subRegion, &interpolationKernel );
        OBLOperatorsKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY >::template launch< POLICY >( subRegion.size(), kernel );
      }
    } );
  }

};

} // namespace oblFluidKernels

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_OBLFLUIDKERNELS_HPP_
