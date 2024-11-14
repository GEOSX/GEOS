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
 * @file SolutionScalingAndCheckingKernelBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONSCALINGANDCHECKINGKERNELBASE_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONSCALINGANDCHECKINGKERNELBASE_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/ElementSubRegionBase.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/**
 * @class SolutionScalingAndCheckingKernelBase
 * @brief Define the kernel for scaling the solution and check its validity
 */
template< typename TYPE >
class SolutionScalingAndCheckingKernelBase
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
  SolutionScalingAndCheckingKernelBase( globalIndex const rankOffset,
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

} // namespace isothermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONSCALINGANDCHECKINGKERNELBASE_HPP
