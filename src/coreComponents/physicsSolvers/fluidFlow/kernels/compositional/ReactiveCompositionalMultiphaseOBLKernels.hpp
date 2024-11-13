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
 * @file ReactiveCompositionalMultiphaseOBLKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_REACTIVECOMPOSITIONALMULTIPHASEOBLKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_REACTIVECOMPOSITIONALMULTIPHASEOBLKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "functions/MultivariableTableFunctionKernels.hpp"
#include "mesh/ElementSubRegionBase.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/ReactiveCompositionalMultiphaseOBLFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"



namespace geos
{

namespace reactiveCompositionalMultiphaseOBLKernels
{

static constexpr real64 minValueForDivision = 1e-10;


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
                      MultivariableTableFunctionStaticKernel< numDofs, numOps > OBLOperatorsTable )
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

    m_OBLOperatorsTable.compute( state, OBLVals, OBLDers );

    // we do not perform derivatives unit conversion here:
    // instead we postpone it till all the derivatives are fully formed, and only then apply the factor only once in 'complete' function
    // scaling the whole system might be even better solution (every pressure column needs to be multiplied by pascalToBarMult)
  }

private:

  // inputs
  MultivariableTableFunctionStaticKernel< numDofs, numOps > m_OBLOperatorsTable;

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
                   MultivariableTableFunction const & function )
  {
    internal::kernelLaunchSelectorEnergySwitch( numPhases, numComponents, enableEnergyBalance, [&] ( auto NP, auto NC, auto E )
    {
      integer constexpr ENABLE_ENERGY = E();
      integer constexpr NUM_PHASES = NP();
      integer constexpr NUM_COMPS = NC();
      integer constexpr NUM_DIMS = ENABLE_ENERGY + NUM_COMPS;
      integer constexpr NUM_OPS  = COMPUTE_NUM_OPS( NUM_PHASES, NUM_COMPS, ENABLE_ENERGY );

      OBLOperatorsKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY >
      kernel( subRegion,
              MultivariableTableFunctionStaticKernel< NUM_DIMS, NUM_OPS >( function.getAxisMinimums(),
                                                                           function.getAxisMaximums(),
                                                                           function.getAxisPoints(),
                                                                           function.getAxisSteps(),
                                                                           function.getAxisStepInvs(),
                                                                           function.getAxisHypercubeMults(),
                                                                           function.getHypercubeData()
                                                                           ) );
      OBLOperatorsKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY >::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }

};

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_PHASES number of phases
 * @tparam NUM_COMPS number of components
 * @tparam ENABLE_ENERGY flag if energy balance equation is assembled
 * @brief Compute accumulation term for an element
 */
template< integer NUM_PHASES, integer NUM_COMPS, bool ENABLE_ENERGY >
class ElementBasedAssemblyKernel
{
public:

  /// Compile time value for the number of phases
  static constexpr integer numPhases = NUM_PHASES;
  /// Compile time value for the number of components
  static constexpr integer numComps = NUM_COMPS;
  /// Compile time value for the energy balance switch
  static constexpr integer enableEnergyBalance = ENABLE_ENERGY;

  /// Compile time value for the number of dimensions
  static constexpr integer numDofs = numComps + enableEnergyBalance;
  // /// Compile time value for the number of operators
  static constexpr integer numOps = COMPUTE_NUM_OPS( NUM_PHASES, NUM_COMPS, ENABLE_ENERGY );

  // component accumulation
  static constexpr integer ACC_OP = 0;
  // kinetic reaction
  static constexpr integer KIN_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases;
  // rock internal energy
  static constexpr integer RE_INTER_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs;

  static constexpr real64 secondsToDaysMult = 1.0 / (60 * 60 * 24);
  static constexpr real64 pascalToBarMult = 1.0 / 1e5;


  /**
   * @brief Constructor
   * @param[in] dt the time step length
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( real64 const dt,
                              globalIndex const rankOffset,
                              string const dofKey,
                              ElementSubRegionBase const & subRegion,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    :
    m_dt( dt * secondsToDaysMult ),
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_referencePoreVolume( subRegion.getField< fields::flow::referencePoreVolume >() ),
    m_referenceRockVolume( subRegion.getField< fields::flow::referenceRockVolume >() ),
    m_rockVolumetricHeatCapacity( subRegion.getField< fields::flow::rockVolumetricHeatCapacity >() ),
    m_rockKineticRateFactor( subRegion.getField< fields::flow::rockKineticRateFactor >() ),
    m_OBLOperatorValues ( subRegion.getField< fields::flow::OBLOperatorValues >()),
    m_OBLOperatorValues_n ( subRegion.getField< fields::flow::OBLOperatorValues_n >()),
    m_OBLOperatorDerivatives ( subRegion.getField< fields::flow::OBLOperatorDerivatives >()),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
    /// Index of the local row corresponding to this element
    localIndex localRow = -1;

    /// Indices of the matrix rows/columns corresponding to the dofs in this element
    globalIndex dofIndices[numDofs]{};

    /// C-array storage for the element local residual vector (all equations except volume balance)
    real64 localResidual[numDofs]{};

    /// C-array storage for the element local Jacobian matrix (all equations except volume balance, all dofs)
    real64 localJacobian[numDofs][numDofs]{};

  };

  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOS_HOST_DEVICE
  inline
  integer elemGhostRank( localIndex const ei ) const
  { return m_elemGhostRank( ei ); }


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {

    // set row index and degrees of freedom indices for this element
    stack.localRow = m_dofNumber[ei] - m_rankOffset;

    for( integer idof = 0; idof < numDofs; ++idof )
    {
      stack.dofIndices[idof] = m_dofNumber[ei] + idof;
    }
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  inline
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack ) const
  {
    arraySlice1d< real64 const, compflow::USD_OBL_VAL - 1 > const & OBLVals = m_OBLOperatorValues[ei];
    arraySlice1d< real64 const, compflow::USD_OBL_VAL - 1 > const & OBLVals_n = m_OBLOperatorValues_n[ei];
    arraySlice2d< real64 const, compflow::USD_OBL_DER - 1 > const & OBLDers = m_OBLOperatorDerivatives[ei];

    // [1] fill diagonal part for both mass (and energy equations if needed, only fluid energy is involved here)
    for( integer c = 0; c < numDofs; ++c )
    {
      stack.localResidual[c] = m_referencePoreVolume[ei] * (OBLVals[ACC_OP + c] - OBLVals_n[ACC_OP + c]);   // acc operators
      // only

      // Add reaction term to diagonal of reservoir cells (here the volume is pore volume or block volume):
      stack.localResidual[c] += (m_referencePoreVolume[ei] + m_referenceRockVolume[ei]) * m_dt * OBLVals[KIN_OP + c] * m_rockKineticRateFactor[ei]; // kinetics

      for( integer v = 0; v < numDofs; ++v )
      {
        stack.localJacobian[c][v] = m_referencePoreVolume[ei] * OBLDers[ACC_OP + c][v];   // der of accumulation term

        // Include derivatives for reaction term if part of reservoir cells:
        stack.localJacobian[c][v] += (m_referencePoreVolume[ei] + m_referenceRockVolume[ei]) * m_dt * OBLDers[KIN_OP + c][v] * m_rockKineticRateFactor[ei];     // derivative
      }
    }

    // + rock energy
    if( enableEnergyBalance )
    {
      stack.localResidual[numDofs-1] += m_referenceRockVolume[ei] * (OBLVals[RE_INTER_OP] - OBLVals_n[RE_INTER_OP]) * m_rockVolumetricHeatCapacity[ei];

      for( integer v = 0; v < numDofs; ++v )
      {
        stack.localJacobian[numDofs-1][v] += m_referenceRockVolume[ei] * OBLDers[RE_INTER_OP][v] * m_rockVolumetricHeatCapacity[ei];
      }   // end of fill offdiagonal part + contribute to diagonal
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  inline
  void complete( localIndex const GEOS_UNUSED_PARAM( ei ),
                 StackVariables & stack ) const
  {
    // add contribution to residual and jacobian into component mass balance equations
    // apply pressure derivative unit conversion
    for( integer i = 0; i < numDofs; ++i )
    {
      stack.localJacobian[i][0] *= pascalToBarMult;
      m_localRhs[stack.localRow + i] += stack.localResidual[i];
      m_localMatrix.addToRow< serialAtomic >( stack.localRow + i,
                                              stack.dofIndices,
                                              stack.localJacobian[i],
                                              numDofs );
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
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

// Time step size
  real64 m_dt;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;

  // views on solid properties
  arrayView1d< real64 const > const m_referencePoreVolume;
  arrayView1d< real64 const > const m_referenceRockVolume;
  arrayView1d< real64 const > const m_rockVolumetricHeatCapacity;
  arrayView1d< real64 const > const m_rockKineticRateFactor;

  // Views on OBL operators and their derivatives
  arrayView2d< real64 const, compflow::USD_OBL_VAL > const m_OBLOperatorValues;
  arrayView2d< real64 const, compflow::USD_OBL_VAL > const m_OBLOperatorValues_n;
  arrayView3d< real64 const, compflow::USD_OBL_DER > const m_OBLOperatorDerivatives;

  // outputs

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

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
   * @param[in] numPhases the number of fluid phases
   * @param[in] numComps the number of fluid components
   * @param[in] enableEnergyBalance flag if energy balance equation is assembled
   * @param[in] dt timestep size
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch(
    integer const numPhases,
    integer const numComps,
    bool const enableEnergyBalance,
    real64 const dt,
    globalIndex const rankOffset,
    string const dofKey,
    ElementSubRegionBase const & subRegion,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs )
  {
    internal::kernelLaunchSelectorEnergySwitch( numPhases, numComps, enableEnergyBalance, [&] ( auto NP, auto NC, auto E )
    {
      integer constexpr ENABLE_ENERGY = E();
      integer constexpr NUM_PHASES = NP();
      integer constexpr NUM_COMPS = NC();

      ElementBasedAssemblyKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY >
      kernel( dt, rankOffset, dofKey, subRegion, localMatrix, localRhs );
      ElementBasedAssemblyKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY >::template launch< POLICY >( subRegion.size(), kernel );
    } );

  }

};

/******************************** FluxComputeKernel ********************************/

/**
 * @brief Base class for FluxComputeKernel that holds all data not dependent
 *        on template parameters (like stencil type and number of components/dofs).
 *
 *        FluxComputeKernel is used for flux terms calculation.
 *        In case mesh geometry/configuration is not changing during simulation,
 *        all connections can be pre-computed, sorted by element, and stored.
 *        Then, ElementBasedAssemblyKernel can be used for flux calculation: every flux will be computed twice,
 *        but data races during matrix and RHS data writes will be avoided, atomic operations are then no longer needed,
 *        and therefore overall performance can significantly improve
 *
 */
class FluxComputeKernelBase
{
public:

  static constexpr real64 pascalToBarMult = 1.0 / 1e5;

  static constexpr real64 secondsToDaysMult = 1.0 / (60 * 60 * 24);

  // transmissibility in DARTS is the same as in Eclipse (Metric):
  // T = c * (k * A) / d, where c is Darcy constant, k is permeability [mD], A is area [m2] and d is distance [m]
  // Darcy constant takes care of unit translation (from SI to Metric), it includes conversion of [s]->[day], [cp->Pa * s], [Pa]->[bar] and
  // [mD->m2]:
  // c = 1e3 [cp]/[Pa*s] * (9.869233e-16) [m2]/[mD] * 86400 [day]/[s] * 1e5 [bar]/[Pa] = 0.008527017312
  // For already forgotten reasons, that constant in DARTS is taken as c[DARTS] = 0.00852671467191601
  // In GEOSX, there is no need in such conversion, as all units are already SI.
  // Therefore, to transform transmissibility in GEOSX to one expected in DARTS,
  // we need to multiply it by the same Darcy constant, but without permeability translation factor:
  // T[DARTS] = T[GEOSX] * c[DARTS] / 9.869233e-16 = T[GEOSX] * 8639693349945.239

  static constexpr real64 transUnitMult = 8639693349945.239;
  static constexpr real64 transDUnitMult = 0.00852671467191601;

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using DofNumberAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  using CompFlowAccessors =
    StencilAccessors< fields::ghostRank,
                      fields::flow::gravityCoefficient,
                      fields::flow::pressure,
                      fields::flow::referencePorosity,
                      fields::flow::rockThermalConductivity,
                      fields::flow::OBLOperatorValues,
                      fields::flow::OBLOperatorDerivatives >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< constitutive::PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofNumberAccessor accessor for data associated with degrees of freedom
   * @param[in] compFlowAccessors accessors for data associated with compositional flow
   * @param[in] permeabilityAccessors accessors for data associated with permeability
   * @param[in] dt time step size
   * @param[in] transMultExp exponent of transmissibility multiplier
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FluxComputeKernelBase( globalIndex const rankOffset,
                         DofNumberAccessor const & dofNumberAccessor,
                         CompFlowAccessors const & compFlowAccessors,
                         PermeabilityAccessors const & permeabilityAccessors,
                         real64 const & dt,
                         real64 const & transMultExp,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs )
    : m_rankOffset( rankOffset ),
    m_dt( dt * secondsToDaysMult ),
    m_transMultExp ( transMultExp ),
    m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
    m_permeability( permeabilityAccessors.get( fields::permeability::permeability {} ) ),
    m_dPerm_dPres( permeabilityAccessors.get( fields::permeability::dPerm_dPressure {} ) ),
    m_referencePorosity( compFlowAccessors.get( fields::flow::referencePorosity {} ) ),
    m_rockThermalConductivity( compFlowAccessors.get( fields::flow::rockThermalConductivity {} ) ),
    m_ghostRank( compFlowAccessors.get( fields::ghostRank {} ) ),
    m_gravCoef( compFlowAccessors.get( fields::flow::gravityCoefficient {} ) ),
    m_pres( compFlowAccessors.get( fields::flow::pressure {} ) ),
    m_OBLOperatorValues ( compFlowAccessors.get( fields::flow::OBLOperatorValues {} ) ),
    m_OBLOperatorDerivatives ( compFlowAccessors.get( fields::flow::OBLOperatorDerivatives {} ) ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

protected:

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Time step size
  real64 const m_dt;

  /// trans multiplier exponent
  real64 m_transMultExp;

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;

  /// Views on permeability
  ElementViewConst< arrayView3d< real64 const > > m_permeability;
  ElementViewConst< arrayView3d< real64 const > > m_dPerm_dPres;
  ElementViewConst< arrayView1d< real64 const > > m_referencePorosity;
  ElementViewConst< arrayView1d< real64 const > > m_rockThermalConductivity;


  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;

  // Primary variables

  /// Views on pressure
  ElementViewConst< arrayView1d< real64 const > > const m_pres;

  // Views on OBL operators and their derivatives
  ElementViewConst< arrayView2d< real64 const, compflow::USD_OBL_VAL > > const m_OBLOperatorValues;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_OBL_DER > > const m_OBLOperatorDerivatives;

  // Residual and jacobian

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;
};

/**
 * @class FluxComputeKernel
 * @tparam NUM_PHASES number of phases
 * @tparam NUM_COMPS number of components
 * @tparam ENABLE_ENERGY flag if energy balance equation is assembled
 * @tparam STENCILWRAPPER wrapper for element stencils
 * @brief Compute flux term for an element
 */
template< integer NUM_PHASES, integer NUM_COMPS, bool ENABLE_ENERGY, typename STENCILWRAPPER >
class FluxComputeKernel : public FluxComputeKernelBase
{
public:

  /// Compile time value for the number of phases
  static constexpr integer numPhases = NUM_PHASES;
  /// Compile time value for the number of components
  static constexpr integer numComps = NUM_COMPS;
  /// Compile time value for the energy balance switch
  static constexpr integer enableEnergyBalance = ENABLE_ENERGY;

  /// Compile time value for the number of primary variables
  static constexpr integer numDofs = numComps + enableEnergyBalance;

  /// Compile time value for the number of equations
  static constexpr integer numEqns = numDofs;

  // /// Compile time value for the number of operators
  static constexpr integer numOps = COMPUTE_NUM_OPS( NUM_PHASES, NUM_COMPS, ENABLE_ENERGY );

  // component flux
  static constexpr integer FLUX_OP = numDofs;
  // diffusion
  static constexpr integer UPSAT_OP = numDofs + numDofs * numPhases;
  static constexpr integer GRAD_OP = numDofs + numDofs * numPhases + numPhases;

  // temperature
  static constexpr integer RE_TEMP_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 1;
  // rock conduction
  static constexpr integer ROCK_COND_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 2;
  // gravity (density)
  static constexpr integer GRAV_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 3;
  // capillary pressure
  static constexpr integer PC_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 3 + numPhases;
  // porosity
  static constexpr integer PORO_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 3 + 2 * numPhases;

  /// Maximum number of elements at the face
  static constexpr localIndex maxNumElems = STENCILWRAPPER::maxNumPointsInFlux;

  /// Maximum number of connections at the face
  static constexpr localIndex maxNumConns = STENCILWRAPPER::maxNumConnections;

  /// Maximum number of points in the stencil
  static constexpr localIndex maxStencilSize = STENCILWRAPPER::maxStencilSize;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor accessor for data associated with degrees of freedom
   * @param[in] compFlowAccessors accessors for data associated with compositional flow
   * @param[in] permeabilityAccessors accessors for data associated with permeability
   * @param[in] dt time step size
   * @param[in] transMultExp exponent of transmissibility multiplier
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FluxComputeKernel( globalIndex const rankOffset,
                     STENCILWRAPPER const & stencilWrapper,
                     DofNumberAccessor const & dofNumberAccessor,
                     CompFlowAccessors const & compFlowAccessors,
                     PermeabilityAccessors const & permeabilityAccessors,
                     real64 const & dt,
                     real64 const & transMultExp,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs )
    : FluxComputeKernelBase( rankOffset,
                             dofNumberAccessor,
                             compFlowAccessors,
                             permeabilityAccessors,
                             dt,
                             transMultExp,
                             localMatrix,
                             localRhs ),
    m_stencilWrapper( stencilWrapper ),
    m_seri( stencilWrapper.getElementRegionIndices() ),
    m_sesri( stencilWrapper.getElementSubRegionIndices() ),
    m_sei( stencilWrapper.getElementIndices() )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {

    /**
     * @brief Constructor for the stack variables
     */
    GEOS_HOST_DEVICE
    StackVariables()
      : dofColIndices( maxNumElems * numDofs ),
      localFlux( numEqns ),
      localFluxJacobian( numEqns, maxNumElems * numDofs )
    {}

    /// Transmissibility
    real64 transmissibility[maxNumConns][2]{};
    /// Thermal Transmissibility
    real64 diffusiveTransmissibility[maxNumConns][2]{};
    /// Derivatives of transmissibility with respect to pressure
    real64 dTrans_dPres[maxNumConns][2]{};

    // Local degrees of freedom and local residual/jacobian

    /// Indices of the matrix rows/columns corresponding to the dofs in this face
    stackArray1d< globalIndex, maxNumElems * numDofs > dofColIndices;

    /// Storage for the face local residual vector
    stackArray1d< real64, numEqns > localFlux;
    /// Storage for the face local Jacobian matrix
    stackArray2d< real64, maxNumElems * numEqns *  numDofs > localFluxJacobian;

  };



  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    // The kernel is only designed for TPFA
    GEOS_ASSERT_EQ( maxNumElems, 2 );
    GEOS_ASSERT_EQ( maxNumConns, 1 );
    GEOS_ASSERT_EQ( maxStencilSize, 2 );

    // set degrees of freedom indices for this face
    for( integer i = 0; i < maxStencilSize; ++i )
    {
      globalIndex const offset = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];

      for( integer jdof = 0; jdof < numDofs; ++jdof )
      {
        stack.dofColIndices[i * numDofs + jdof] = offset + jdof;
      }
    }
  }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  inline
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {



    // first, compute the transmissibilities at this face
    // though weights are recomputed every time, the perm derivative is not used:
    // perm/poro changes are treated differently in OBL: via m_transMultExp and PORO operator

    // in case mesh geometry/configuration is not changing during simulation, weights can be computed as a preprocessing step
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );

    m_stencilWrapper.computeWeights( iconn,
                                     stack.diffusiveTransmissibility,
                                     stack.dTrans_dPres );


    // As a first iteration, do everything in TPFA style: i is the left (minus) element, j is the right(plus) element.
    // Consider inflow and outflow to and from element i.

    localIndex const erI  = m_seri( iconn, 0 );
    localIndex const esrI = m_sesri( iconn, 0 );
    localIndex const eiI  = m_sei( iconn, 0 );

    localIndex const erJ  = m_seri( iconn, 1 );
    localIndex const esrJ = m_sesri( iconn, 1 );
    localIndex const eiJ  = m_sei( iconn, 1 );

    arraySlice1d< real64 const, compflow::USD_OBL_VAL - 1 > const & OBLValsI = m_OBLOperatorValues[erI][esrI][eiI];
    arraySlice2d< real64 const, compflow::USD_OBL_DER - 1 > const & OBLDersI = m_OBLOperatorDerivatives[erI][esrI][eiI];

    arraySlice1d< real64 const, compflow::USD_OBL_VAL - 1 > const & OBLValsJ = m_OBLOperatorValues[erJ][esrJ][eiJ];
    arraySlice2d< real64 const, compflow::USD_OBL_DER - 1 > const & OBLDersJ = m_OBLOperatorDerivatives[erJ][esrJ][eiJ];

    // single gravity coefficient and transimissibility values correspond to the current connection
    // we need to apply pascal -> bar conversion to gravity coefficient, the same way as in DARTS
    real64 const gravCoef = (m_gravCoef[erI][esrI][eiI] - m_gravCoef[erJ][esrJ][eiJ]) * pascalToBarMult;
    real64 const trans = stack.transmissibility[0][0] * transUnitMult;
    real64 const transD = stack.diffusiveTransmissibility[0][0] * transDUnitMult;


    real64 transMult = 1;
    real64 transMultD = 0;
    if( m_transMultExp > 0 )
    {
      // Calculate transmissibility multiplier:

      // Take average interface porosity:
      real64 const poroAverage = (OBLValsI[PORO_OP] + OBLValsJ[PORO_OP]) * 0.5;
      transMultD = m_transMultExp * pow( poroAverage, m_transMultExp - 1 ) * 0.5;
      transMult = pow( poroAverage, m_transMultExp );
    }


    // apply [Pa]->[bar] conversion
    real64 const pDiff = (m_pres[erJ][esrJ][eiJ] - m_pres[erI][esrI][eiI]) * pascalToBarMult;


    // [2] fill offdiagonal part + contribute to diagonal, only fluid part is considered in energy equation
    for( integer p = 0; p < numPhases; ++p )
    {
      // loop over number of phases for convective operator

      // calculate gravity term for phase p
      real64 const averageDensity = (OBLValsI[GRAV_OP + p] + OBLValsJ[GRAV_OP + p]) / 2;

      // p = 1 means oil phase, it's reference phase. pw=po-pcow, pg=po-(-pcog).

      real64 const phasePDiff = pDiff + averageDensity * gravCoef - OBLValsJ[PC_OP + p] + OBLValsI[PC_OP + p];

      // note: pressure conversion factor is only applied for RHS and not for Jacobian:
      // instead we postpone it till all the derivatives are fully formed, and only then apply the factor only once in 'complete' function
      // scaling the whole system might be even better solution (every pressure column needs to be multiplied by pascalToBarMult)

      // calculate partial derivatives for gravity and capillary terms
      real64 gravPcDerI[numDofs];
      real64 gravPcDerJ[numDofs];
      for( integer v = 0; v < numDofs; ++v )
      {
        gravPcDerI[v] = -(OBLDersI[GRAV_OP + p][v]) * gravCoef / 2 - OBLDersI[PC_OP + p][v];
        gravPcDerJ[v] = -(OBLDersJ[GRAV_OP + p][v]) * gravCoef / 2 + OBLDersJ[PC_OP + p][v];
      }

      real64 phaseGammaPDiff = transMult * trans * m_dt * phasePDiff;

      arraySlice1d< real64 const, compflow::USD_OBL_VAL - 1 > const & OBLValsUp = (phasePDiff < 0) ? OBLValsI : OBLValsJ;
      arraySlice2d< real64 const, compflow::USD_OBL_DER - 1 > const & OBLDersUp = (phasePDiff < 0) ? OBLDersI : OBLDersJ;
      integer const upOffset = (phasePDiff < 0) ? 0 : numDofs;

      // mass and energy outflow with effect of gravity and capillarity
      for( integer c = 0; c < numDofs; ++c )
      {
        real64 compFlux = transMult * trans * m_dt * OBLValsUp[FLUX_OP + p * numDofs + c];

        stack.localFlux[c] -= phasePDiff * compFlux;       // flux operators only
        for( integer v = 0; v < numDofs; ++v )
        {
          stack.localFluxJacobian[c][v + upOffset] -= (phaseGammaPDiff * OBLDersUp[FLUX_OP + p * numDofs + c][v] +
                                                       trans * m_dt * phasePDiff * transMultD * OBLDersUp[PORO_OP][v] * OBLValsUp[FLUX_OP + p * numDofs + c]);
          stack.localFluxJacobian[c][v] += compFlux * gravPcDerI[v];
          stack.localFluxJacobian[c][numDofs + v] += compFlux * gravPcDerJ[v];

          // contribution from pressure derivative related to pDiff in phasePDiff * compFlux term
          if( v == 0 )
          {
            stack.localFluxJacobian[c][v] += compFlux;
            stack.localFluxJacobian[c][numDofs + v] -= compFlux;
          }
        }
      }

    }     // end of loop over number of phases for convective operator with gravity and capillarity

    // [3] Additional diffusion code here:   (phi_p * S_p) * (rho_p * D_cp * Delta_x_cp)  or (phi_p * S_p) * (kappa_p * Delta_T)

    real64 const poroAverage = (m_referencePorosity[erI][esrI][eiI] + m_referencePorosity[erJ][esrJ][eiJ]) * 0.5;     // diffusion term
                                                                                                                      // depends on total
    // porosity!


    // Add diffusion term to the residual:
    for( integer c = 0; c < numDofs; ++c )
    {
      for( integer p = 0; p < numPhases; ++p )
      {
        real64 phaseCompGrad = OBLValsJ[GRAD_OP + c * numPhases + p] - OBLValsI[GRAD_OP + c * numPhases + p];

        arraySlice1d< real64 const, compflow::USD_OBL_VAL - 1 > const & OBLValsUp = (phaseCompGrad < 0) ? OBLValsI : OBLValsJ;
        arraySlice2d< real64 const, compflow::USD_OBL_DER - 1 > const & OBLDersUp = (phaseCompGrad < 0) ? OBLDersI : OBLDersJ;
        integer const upOffset = (phaseCompGrad < 0) ? 0 : numDofs;

        // use upstream quantity from cell i for compressibility and saturation (mass or energy):
        real64 phaseGammaDDiff = m_dt * transD * poroAverage * OBLValsUp[UPSAT_OP + p];

        stack.localFlux[c] -= phaseGammaDDiff * phaseCompGrad;         // diffusion term

        // Add diffusion terms to Jacobian:
        for( integer v = 0; v < numDofs; ++v )
        {
          stack.localFluxJacobian[c][v] += phaseGammaDDiff * OBLDersI[GRAD_OP + c * numPhases + p][v];
          stack.localFluxJacobian[c][numDofs + v] -= phaseGammaDDiff * OBLDersJ[GRAD_OP + c * numPhases + p][v];

          stack.localFluxJacobian[c][v + upOffset] -= phaseCompGrad * m_dt * transD * poroAverage * OBLDersUp[UPSAT_OP + p][v];
        }

      }
    }

    // [4] add rock conduction
    if( enableEnergyBalance )
    {
      real64 const tDiff = OBLValsJ[RE_TEMP_OP] - OBLValsI[RE_TEMP_OP];
      real64 const gammaTDiff = transD * m_dt * tDiff;

      arraySlice1d< real64 const, compflow::USD_OBL_VAL - 1 > const & OBLValsUp = (tDiff < 0) ? OBLValsI : OBLValsJ;
      arraySlice2d< real64 const, compflow::USD_OBL_DER - 1 > const & OBLDersUp = (tDiff < 0) ? OBLDersI : OBLDersJ;
      integer const upOffset = (tDiff < 0) ? 0 : numDofs;
      localIndex const erUp = (tDiff < 0) ? erI : erJ;
      localIndex const esrUp = (tDiff < 0) ? esrI : esrJ;
      localIndex const eiUp = (tDiff < 0) ? eiI : eiJ;

      real64 const poroUp = m_referencePorosity[erUp][esrUp][eiUp];
      real64 const rockCondUp = m_rockThermalConductivity[erUp][esrUp][eiUp];

      // rock heat transfers flows from cell i to j
      stack.localFlux[numComps] -= gammaTDiff * OBLValsUp[ROCK_COND_OP] * (1 - poroUp) * rockCondUp;
      for( integer v = 0; v < numDofs; ++v )
      {
        stack.localFluxJacobian[numComps][v + upOffset] -= gammaTDiff * OBLDersUp[ROCK_COND_OP][v] * (1 - poroUp) * rockCondUp;
        // the last variable - T
        if( v == numDofs - 1 )
        {
          real64 const tFlux = transD * m_dt * OBLValsUp[ROCK_COND_OP] * (1 - poroUp) * rockCondUp;
          stack.localFluxJacobian[numComps][v] += tFlux;
          stack.localFluxJacobian[numComps][numDofs + v] -= tFlux;
        }
      }
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  inline
  void complete( localIndex const iconn,
                 StackVariables & stack ) const
  {
    // Add to residual/jacobian
    for( integer i = 0; i < maxNumElems; ++i )
    {
      // during first loop iteration (i == 0), we fill the equations for element i as is,
      // only scaling the pressure derivative to account for unit conversion
      if( i == 0 )
      {
        for( integer id = 0; id < numDofs; ++id )
        {
          stack.localFluxJacobian[id][0] *= pascalToBarMult;
          stack.localFluxJacobian[id][numDofs] *= pascalToBarMult;
        }
      }
      // during second (and the last - TPFA) loop iteration (i ==1), we fill now the equations for element j.
      // here, we need to multiply all values by -1, since the flow direction changes w.r.t element
      else
      {
        for( integer id = 0; id < numDofs; ++id )
        {
          stack.localFlux[id] *= -1;
          for( integer v = 0; v < 2 * numDofs; ++v )
          {
            stack.localFluxJacobian[id][v] *= -1;
          }
        }
      }

      if( m_ghostRank[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
        GEOS_ASSERT_GE( localRow, 0 );
        GEOS_ASSERT_GE( m_localMatrix.numRows(), localRow + numEqns );

        for( integer id = 0; id < numDofs; ++id )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow + id], stack.localFlux[id] );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
            ( localRow + id,
            stack.dofColIndices.data(),
            stack.localFluxJacobian[id].dataIfContiguous(),
            maxNumElems * numDofs );
        }
      }
    }
  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numConnections the number of connections
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numConnections,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    if( numConnections )
    {
      forAll< POLICY >( numConnections, [=] GEOS_HOST_DEVICE ( localIndex const iconn )
      {
        typename KERNEL_TYPE::StackVariables stack;

        kernelComponent.setup( iconn, stack );
        kernelComponent.computeFlux( iconn, stack );
        kernelComponent.complete( iconn, stack );
      } );
    }
  }

private:

  // Stencil information

  /// Reference to the stencil wrapper
  STENCILWRAPPER const m_stencilWrapper;

  /// Connection to element maps
  typename STENCILWRAPPER::IndexContainerViewConstType const m_seri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sesri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sei;
};

/**
 * @class FluxComputeKernelFactory
 */
class FluxComputeKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam STENCILWRAPPER the type of the stencil wrapper
   * @param[in] numPhases the number of fluid phases
   * @param[in] numComps the number of fluid components
   * @param[in] enableEnergyBalance flag if energy balance equation is assembled
   * @param[in] transMultExp the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numPhases,
                   integer const numComps,
                   bool const enableEnergyBalance,
                   real64 const & transMultExp,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    internal::kernelLaunchSelectorEnergySwitch( numPhases, numComps, enableEnergyBalance, [&] ( auto NP, auto NC, auto E )
    {
      integer constexpr ENABLE_ENERGY = E();
      integer constexpr NUM_PHASES = NP();
      integer constexpr NUM_COMPS = NC();

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      using KERNEL_TYPE = FluxComputeKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY, STENCILWRAPPER >;
      typename KERNEL_TYPE::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

      KERNEL_TYPE kernel( rankOffset, stencilWrapper, dofNumberAccessor, compFlowAccessors, permeabilityAccessors,
                          dt, transMultExp, localMatrix, localRhs );
      KERNEL_TYPE::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};


/******************************** ResidualNormKernels ********************************/

struct ResidualNormKernel
{

  template< typename POLICY >
  static void launch( arrayView1d< real64 const > const & localResidual,
                      globalIndex const rankOffset,
                      integer const numDofs,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< integer const > const & ghostRank,
                      arrayView1d< real64 const > const & refPoreVolume,
                      arrayView2d< real64 const, compflow::USD_OBL_VAL > const & OBLOperatorValues_n,
                      real64 & localResidualNorm )
  {
    RAJA::ReduceSum< ReducePolicy< POLICY >, real64 > localSum( 0.0 );

    forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        localIndex const localRow = dofNumber[ei] - rankOffset;

        for( integer idof = 0; idof < numDofs; ++idof )
        {
          real64 const val = localResidual[localRow + idof] / (OBLOperatorValues_n[ei][idof] * refPoreVolume[ei]);
          localSum += val * val;
        }
      }
    } );
    localResidualNorm += localSum.get();
  }
};

struct ResidualDARTSL2NormKernel
{

  template< typename POLICY >
  static void launch( arrayView1d< real64 const > const & localResidual,
                      globalIndex const rankOffset,
                      integer const numDofs,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< integer const > const & ghostRank,
                      arrayView1d< real64 const > const & refPoreVolume,
                      arrayView2d< real64 const, compflow::USD_OBL_VAL > const & OBLOperatorValues,
                      real64 & localResidualNorm )
  {
    real64 constexpr eps = minValueForDivision;

    for( integer idof = 0; idof < numDofs; ++idof )
    {
      RAJA::ReduceSum< ReducePolicy< POLICY >, real64 > localResSum( 0.0 ), localNormSum( 0.0 );

      forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          localIndex const localRow = dofNumber[ei] - rankOffset;

          localResSum += localResidual[localRow + idof] * localResidual[localRow + idof];
          localNormSum +=  OBLOperatorValues[ei][idof] * OBLOperatorValues[ei][idof] * refPoreVolume[ei] * refPoreVolume[ei];
        }
      } );

      real64 const normalizer = ( localNormSum.get() > eps ) ? localNormSum.get() : eps;
      real64 const norm = sqrt( localResSum.get() / normalizer );
      if( localResidualNorm < norm )
      {
        localResidualNorm = norm;
      }
    }
  }
};


/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY >
  static localIndex
  launch( arrayView1d< real64 const > const & localSolution,
          globalIndex const rankOffset,
          localIndex const numComps,
          bool const enableThermalBalance,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & pres,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac,
          arrayView1d< real64 const > const & temp,
          integer const allowOBLChopping,
          real64 const scalingFactor )
  {
    real64 constexpr eps = minValueForDivision;

    RAJA::ReduceMin< ReducePolicy< POLICY >, integer > check( 1 );

    forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        localIndex const localRow = dofNumber[ei] - rankOffset;

        real64 const newPres = pres[ei] + scalingFactor * localSolution[localRow];
        check.min( newPres >= 0.0 );

        // if OBL chopping is not allowed, the time step fails if a component fraction is negative
        // otherwise, we just check that the total fraction is positive, and out-of-OBL-bounds component fractions
        // will be chopped (i.e., set to minimum OBL limit) in ApplySystemSolution)
        if( !allowOBLChopping )
        {
          for( integer ic = 0; ic < numComps - 1; ++ic )
          {
            real64 const newCompFrac = compFrac[ei][ic] + scalingFactor * localSolution[localRow + ic + 1];
            check.min( newCompFrac >= 0.0 );
          }
        }
        else
        {
          real64 fracSum = 0.0;
          for( integer ic = 0; ic < numComps - 1; ++ic )
          {
            real64 const newCompFrac = compFrac[ei][ic] + scalingFactor * localSolution[localRow + ic + 1];
            fracSum += (newCompFrac > 0.0) ? newCompFrac : 0.0;
          }
          check.min( fracSum >= eps );
        }

        if( enableThermalBalance )
        {
          real64 const newTemp = temp[ei] + scalingFactor * localSolution[localRow + numComps];
          check.min( newTemp >= 0.0 );
        }

      }
    } );
    return check.get();
  }

};

} // namespace ReactiveCompositionalMultiphaseOBLKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_REACTIVECOMPOSITIONALMULTIPHASEOBLKERNELS_HPP
