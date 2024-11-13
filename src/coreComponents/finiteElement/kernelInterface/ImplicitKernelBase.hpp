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

#include "KernelBase.hpp"

/**
 * @file ImplicitKernelBase.hpp
 */

#ifndef GEOS_FINITEELEMENT_IMPLICITKERNELBASE_HPP_
#define GEOS_FINITEELEMENT_IMPLICITKERNELBASE_HPP_



namespace geos
{

namespace finiteElement
{

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @class ImplicitKernelBase
 * @brief Define the base interface for implicit finite element kernels.
 * @copydoc geos::finiteElement::KernelBase
 *
 * ### ImplicitKernelBase Description
 * Provides a common base for kernels that require the assembly of a system of
 * equations. The types required to assemble the system, such as DOF
 * information, the Matrix and Vector object, etc., are declared and set here.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class ImplicitKernelBase : public KernelBase< SUBREGION_TYPE,
                                              CONSTITUTIVE_TYPE,
                                              FE_TYPE,
                                              NUM_DOF_PER_TEST_SP,
                                              NUM_DOF_PER_TRIAL_SP >
{
public:
  /// Alias for the base class. (i.e. #geos::finiteElement::KernelBase)
  using Base = KernelBase< SUBREGION_TYPE,
                           CONSTITUTIVE_TYPE,
                           FE_TYPE,
                           NUM_DOF_PER_TEST_SP,
                           NUM_DOF_PER_TRIAL_SP >;

  using Base::maxNumTestSupportPointsPerElem;
  using Base::maxNumTrialSupportPointsPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_finiteElementSpace;


  /**
   * @brief Constructor
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param inputDofNumber The dof number for the primary field.
   * @param rankOffset dof index offset of current rank
   * @param inputMatrix Reference to the Jacobian matrix.
   * @param inputRhs Reference to the RHS vector.
   * @param inputDt The timestep for the physics update.
   * @copydoc geos::finiteElement::KernelBase::KernelBase
   */
  ImplicitKernelBase( NodeManager const & nodeManager,
                      EdgeManager const & edgeManager,
                      FaceManager const & faceManager,
                      localIndex const targetRegionIndex,
                      SUBREGION_TYPE const & elementSubRegion,
                      FE_TYPE const & finiteElementSpace,
                      CONSTITUTIVE_TYPE & inputConstitutiveType,
                      arrayView1d< globalIndex const > const & inputDofNumber,
                      globalIndex const rankOffset,
                      CRSMatrixView< real64, globalIndex const > const & inputMatrix,
                      arrayView1d< real64 > const & inputRhs,
                      real64 const inputDt ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_dofNumber( inputDofNumber ),
    m_dofRankOffset( rankOffset ),
    m_matrix( inputMatrix ),
    m_rhs( inputRhs ),
    m_dt( inputDt )
  {
    FiniteElementBase::initialize< FE_TYPE >( nodeManager,
                                              edgeManager,
                                              faceManager,
                                              elementSubRegion,
                                              m_meshData );
    GEOS_UNUSED_VAR( targetRegionIndex );
  }


  //***************************************************************************
  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
    /// The number of rows in the pre-allocated element local jacobian matrix (upper bound for numRows).
    static constexpr int maxNumRows = maxNumTestSupportPointsPerElem *numDofPerTestSupportPoint;

    /// The number of columns in the pre-allocated element local jacobian matrix (upper bound for numCols).
    static constexpr int maxNumCols = maxNumTrialSupportPointsPerElem *numDofPerTrialSupportPoint;

    /**
     * Default constructor
     */
    GEOS_HOST_DEVICE
    StackVariables():
      localRowDofIndex{ 0 },
      localColDofIndex{ 0 },
      localResidual{ 0.0 },
      localJacobian{ {0.0} }
    {
      for( int ii = 0; ii < maxNumRows; ++ii )
      {
        for( int jj = 0; jj < maxNumCols; ++jj )
        {
          localJacobian[ii][jj] = 0.0;
        }
      }
    }

    /// The actual number of rows in the element local jacobian matrix (<= maxNumRows).
    localIndex numRows;

    /// The actual number of columns in the element local jacobian matrix (<= maxNumCols).
    localIndex numCols;

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localRowDofIndex[maxNumRows];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex localColDofIndex[maxNumCols];

    /// C-array storage for the element local residual vector.
    real64 localResidual[maxNumRows];

    /// C-array storage for the element local Jacobian matrix.
    real64 localJacobian[maxNumRows][maxNumCols];

    /// Stack variables needed for the underlying FEM type
    typename FE_TYPE::StackVariables feStack;
  };
  //***************************************************************************

  /**
   * @copydoc geos::finiteElement::KernelBase::setup
   *
   * ### ImplicitKernelBase::setup() Description
   *
   * In this implementation, the element local Row and Column DOF stack arrays
   * are filled for when we fill the global matrix and rhs.
   *
   * @note This seems like a waste of register space. We should do this in
   *       complete() unless we actually need these dof somewhere else in the kernel.
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    m_finiteElementSpace.template setup< FE_TYPE >( k, m_meshData, stack.feStack );
    localIndex numTestSupportPoints = m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
    localIndex numTrialSupportPoints = numTestSupportPoints;
    stack.numRows = numTestSupportPoints * numDofPerTestSupportPoint;
    stack.numCols = numTrialSupportPoints * numDofPerTrialSupportPoint;
    for( localIndex a=0; a<numTestSupportPoints; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes[k][a];
      for( int i=0; i<numDofPerTestSupportPoint; ++i )
      {
        stack.localRowDofIndex[a*numDofPerTestSupportPoint+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

    for( localIndex a=0; a<numTrialSupportPoints; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes[k][a];
      for( int i=0; i<numDofPerTrialSupportPoint; ++i )
      {
        stack.localColDofIndex[a*numDofPerTrialSupportPoint+i] = m_dofNumber[localNodeIndex]+i;
      }
    }
  }


protected:
  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_dofNumber;

  /// The global rank offset
  globalIndex const m_dofRankOffset;

  /// The global Jacobian matrix.
  CRSMatrixView< real64, globalIndex const > const m_matrix;

  /// The global residaul vector.
  arrayView1d< real64 > const m_rhs;

  /// Data structure containing mesh data used to setup the finite element
  typename FE_TYPE::template MeshData< SUBREGION_TYPE > m_meshData;

  /// time increment
  real64 const m_dt; ///TODO: Consider moving to finite element kernel base?

};

}
}



#endif /* GEOS_FINITEELEMENT_IMPLICITKERNELBASE_HPP_ */
