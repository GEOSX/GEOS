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
 * @file AccumulationZFormulationKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_ACCUMULATIONZFORMULATIONKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_ACCUMULATIONZFORMULATIONKERNEL_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "mesh/ElementSubRegionBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/KernelLaunchSelectors.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** AccumulationKernel ********************************/

/**
 * @class AccumulationZFormulationKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of accumulation and volume balance
 */
template< integer NUM_COMP, integer NUM_DOF>
class AccumulationZFormulationKernel
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
  AccumulationZFormulationKernel( localIndex const numPhases,
                      globalIndex const rankOffset,
                      string const dofKey,
                      ElementSubRegionBase const & subRegion,
                      constitutive::MultiFluidBase const & fluid,
                      constitutive::CoupledSolidBase const & solid,
                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                      arrayView1d< real64 > const & localRhs,
                      BitFlags< KernelFlags > const KernelFlags )
    : m_numPhases( numPhases ),
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_volume( subRegion.getElementVolume() ),
    m_porosity( solid.getPorosity() ),
    m_dPoro_dPres( solid.getDporosity_dPressure() ),
    m_phaseVolFrac( subRegion.getField< fields::flow::phaseVolumeFraction >() ),
    m_dPhaseVolFrac( subRegion.getField< fields::flow::dPhaseVolumeFraction >() ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() ),
    m_phaseCompFrac( fluid.phaseCompFraction() ),
    m_dPhaseCompFrac( fluid.dPhaseCompFraction() ),
    m_phaseFrac( fluid.phaseFraction() ),
    m_dPhaseFrac( fluid.dPhaseFraction() ),
    m_compDens( subRegion.getField< fields::flow::globalCompDensity >() ),
    m_compFrac( subRegion.getField< fields::flow::globalCompFraction >() ),
    m_compAmount_n( subRegion.getField< fields::flow::compAmount_n >() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs ),
    m_KernelFlags( KernelFlags )
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
    using Deriv = constitutive::multifluid::DerivativeOffset;

    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const compDens = m_compDens[ei];
    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const compFrac = m_compFrac[ei];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseDens = m_phaseDens[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseDens = m_dPhaseDens[ei][0];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseFrac = m_phaseFrac[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseFrac = m_dPhaseFrac[ei][0];

    // calculate total component density from volume balance
    // TO DO: create separate kernel
    real64 totalDensity = 0.0;
    real64 dtotalDensity_dP = 0.0;
    real64 dtotalDensity_dC[numComp]{};

    for( integer ic = 0; ic < numComp; ++ic )
    {
      totalDensity += compDens[ic];
      dtotalDensity_dC[ic] = 0.0;
    }

    //real64 const sum_phaseFrac_phaseDens = 0.0;
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      real64 const phaseDensInv = 1.0 / phaseDens[ip];
      //sum_phaseFrac_phaseDens += phaseFrac[ip] * phaseDensInv;

      dtotalDensity_dP += phaseDensInv * totalDensity * totalDensity *
        (phaseFrac[ip] * phaseDensInv * dPhaseDens[ip][Deriv::dP] - dPhaseFrac[ip][Deriv::dP]);

      for( integer jc = 0; jc < numComp; ++jc )
      {
        dtotalDensity_dC[jc] += phaseDensInv * totalDensity * totalDensity *
          (phaseFrac[ip] * phaseDensInv * dPhaseDens[ip][Deriv::dC+jc] - dPhaseFrac[ip][Deriv::dC+jc]);
      }
    }

    std::cout << "Cell [" << ei << "]" << std::endl;
    for( integer ic = 0; ic < numComp; ++ic )
      std::cout << "rho[" << ic << "] = " << compDens[ic] << std::endl;
    for( integer ic = 0; ic < numComp; ++ic )
      std::cout << "z[" << ic << "] = " << compFrac[ic] << std::endl;  
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      std::cout << "nu[" << ip << "] = " << phaseFrac[ip] << std::endl; 
      std::cout << "dnu[" << ip << "]_dP = " << dPhaseFrac[ip][Deriv::dP] << std::endl;   
      for( integer ic = 0; ic < numComp; ++ic )
        std::cout << "dnu[" << ip << "]_dz[" << ic << "] = " << dPhaseFrac[ip][Deriv::dC+ic] << std::endl;  
    }

    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      std::cout << "phaseDens[" << ip << "] = " << phaseDens[ip] << std::endl; 
      std::cout << "dphaseDens[" << ip << "]_dP = " << dPhaseDens[ip][Deriv::dP] << std::endl;   
      for( integer ic = 0; ic < numComp; ++ic )
        std::cout << "dphaseDens[" << ip << "]_dz[" << ic << "] = " << dPhaseDens[ip][Deriv::dC+ic] << std::endl;  
    }
      
   

    std::cout << "rhoT = " << totalDensity << std::endl;  
    std::cout << "drhoT_dP = " << dtotalDensity_dP << std::endl;
    for( integer ic = 0; ic < numComp; ++ic )
      std::cout << "drhoT_dC[" << ic << "] = " << dtotalDensity_dC[ic] << std::endl;
    
    // ic - index of component whose conservation equation is assembled
    // (i.e. row number in local matrix)
    for( integer ic = 0; ic < numComp; ++ic )
    {
      real64 const compAmount = stack.poreVolume * totalDensity * compFrac[ic];
      real64 const compAmount_n = m_compAmount_n[ei][ic];

      stack.localResidual[ic] += compAmount - compAmount_n;

      // derivatives with respect to pressure (p)
      real64 const dCompAmount_dP = compFrac[ic] * (stack.dPoreVolume_dPres * totalDensity + stack.poreVolume * dtotalDensity_dP);
      stack.localJacobian[ic][0] += dCompAmount_dP;

      // derivatives with respect to global component fraction (zc)
      for( integer jc = 0; jc < numComp; ++jc )
      {
        real64 dCompAmount_dC;
        if (ic == jc)
          dCompAmount_dC = stack.poreVolume * (totalDensity + dtotalDensity_dC[jc] * compFrac[ic]);
        else
          dCompAmount_dC = stack.poreVolume * (dtotalDensity_dC[jc] * compFrac[ic]);

        stack.localJacobian[ic][jc + 1] += dCompAmount_dC;
      }
    }
  
    for( integer ic = 0; ic < numComp; ++ic )
    {
      std::cout << "Raccum[" << ic << "] = " << stack.localResidual[ic] << std::endl;
      std::cout << "dRaccum[" << ic << "]_dP = " << stack.localJacobian[ic][0] << std::endl;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        std::cout << "dRaccum[" << ic << "]_dz[" << jc << "] = " << stack.localJacobian[ic][jc+1] << std::endl;
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

    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > compFrac = m_compFrac[ei];

    real64 oneMinusCompFracSum = 1.0;

    // sum contributions to component accumulation from each component
    for( integer ic = 0; ic < numComp; ++ic )
    {
      oneMinusCompFracSum -= compFrac[ic];
      // no derivatives w.r.t pressure

      // derivative w.r.t component fractions are unity: dzc_dzc = 1
      stack.localJacobian[numComp][ic+1] -= 1;
    }

    // scale componentFraction-based volume balance by pore volume (for better scaling w.r.t. other equations)
    stack.localResidual[numComp] = stack.poreVolume * oneMinusCompFracSum;
    for( integer idof = 0; idof < numDof; ++idof )
    {
      stack.localJacobian[numComp][idof] *= stack.poreVolume;
    }
    stack.localJacobian[numComp][0] += stack.dPoreVolume_dPres * oneMinusCompFracSum;
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

    if( m_KernelFlags.isSet( KernelFlags::TotalMassEquation ) )
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

  /// Views on phase fractions
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseFrac;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const m_dPhaseFrac;

  /// Views on the phase densities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseDens;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const m_dPhaseDens;

  /// Views on the phase component fraction
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > const m_phaseCompFrac;
  arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > const m_dPhaseCompFrac;

  // View on component densities and component fractions 
  arrayView2d< real64 const, compflow::USD_COMP > m_compDens;
  arrayView2d< real64 const, compflow::USD_COMP > m_compFrac;

  // View on component amount (mass/moles) from previous time step
  arrayView2d< real64 const, compflow::USD_COMP > m_compAmount_n;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  BitFlags< KernelFlags > const m_KernelFlags;
};

/**
 * @class AccumulationZFormulationKernelFactory
 */
class AccumulationZFormulationKernelFactory
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

      BitFlags< KernelFlags > KernelFlags;
      if( useTotalMassEquation )
        KernelFlags.set( KernelFlags::TotalMassEquation );
      if( useSimpleAccumulation )
        KernelFlags.set( KernelFlags::SimpleAccumulation );

      AccumulationZFormulationKernel< NUM_COMP, NUM_DOF> kernel( numPhases, rankOffset, dofKey, subRegion,
                                                      fluid, solid, localMatrix, localRhs, KernelFlags );
      AccumulationZFormulationKernel< NUM_COMP, NUM_DOF>::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }

};

} // namespace isothermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_ACCUMULATIONZFORMULATIONKERNEL_HPP
