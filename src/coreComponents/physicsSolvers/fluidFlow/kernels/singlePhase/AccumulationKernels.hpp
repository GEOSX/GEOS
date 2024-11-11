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
 * @file AccumulationKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_ACCUMULATIONKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_ACCUMULATIONKERNELS_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** AccumulationKernel ********************************/

/**
 * @class AccumulationKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation
 */
template< typename SUBREGION_TYPE, integer NUM_DOF >
class AccumulationKernel
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
  AccumulationKernel( globalIndex const rankOffset,
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
 * @class SurfaceElementAccumulationKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation in SurfaceElementSubRegion
 */
class SurfaceElementAccumulationKernel : public AccumulationKernel< SurfaceElementSubRegion, 1 >
{

public:

  using Base = AccumulationKernel< SurfaceElementSubRegion, 1 >;

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
  SurfaceElementAccumulationKernel( globalIndex const rankOffset,
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
 * @class AccumulationKernelFactory
 */
class AccumulationKernelFactory
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
  template< typename POLICY, typename SUBREGION_TYPE >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const dofKey,
                   SUBREGION_TYPE const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    if constexpr ( std::is_base_of_v< CellElementSubRegion, SUBREGION_TYPE > )
    {
      integer constexpr NUM_DOF = 1;
      AccumulationKernel< CellElementSubRegion, NUM_DOF > kernel( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
      AccumulationKernel< CellElementSubRegion, NUM_DOF >::template launch< POLICY >( subRegion.size(), kernel );
    }
    else if constexpr ( std::is_base_of_v< SurfaceElementSubRegion, SUBREGION_TYPE > )
    {
      SurfaceElementAccumulationKernel kernel( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
      SurfaceElementAccumulationKernel::launch< POLICY >( subRegion.size(), kernel );
    }
    else
    {
      GEOS_UNUSED_VAR( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
      GEOS_ERROR( "Unsupported subregion type: " << typeid(SUBREGION_TYPE).name() );
    }
  }

};

} // namespace singlePhaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_ACCUMULATIONKERNELS_HPP
