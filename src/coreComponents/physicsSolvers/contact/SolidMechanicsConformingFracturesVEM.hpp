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
 * @file SolidMechanicsConformingFracturesVEM.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSCONFORMINGFRACTURESVEM_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSCONFORMINGFRACTURESVEM_HPP_

#include "physicsSolvers/contact/ContactSolverBase.hpp"

namespace geos
{

class SolidMechanicsLagrangianFEM;

class SolidMechanicsConformingFracturesVEM : public ContactSolverBase
{
public:

  SolidMechanicsConformingFracturesVEM( const string & name,
                           Group * const parent );

  ~SolidMechanicsConformingFracturesVEM() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "SolidMechanicsConformingFracturesVEM";
  }

  /// String used to form the solverName used to register single-physics solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "SolidMechanicsConformingFracturesVEM"; }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & MeshBodies ) override final;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = true ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override final;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( real64 const & time,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64
  setNextDt( real64 const & currentDt,
             DomainPartition & domain ) override;

  void updateState( DomainPartition & domain ) override final;

  void assembleForceResidualDerivativeWrtTraction( MeshLevel const & mesh,
                                                   arrayView1d< string const > const & regionNames,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs );

  void assembleTractionResidualDerivativeWrtDisplacementAndTraction( MeshLevel const & mesh,
                                                                     arrayView1d< string const > const & regionNames,
                                                                     DofManager const & dofManager,
                                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                     arrayView1d< real64 > const & localRhs );

  void assembleStabilization( MeshLevel const & mesh,
                              NumericalMethodsManager const & numericalMethodManager,
                              DofManager const & dofManager,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs );

  string const & getContactRelationName() const { return m_contactRelationName; }

  bool resetConfigurationToDefault( DomainPartition & domain ) const override final;

  bool updateConfiguration( DomainPartition & domain ) override final;

  bool isFractureAllInStickCondition( DomainPartition const & domain ) const;

  void computeRotationMatrices( DomainPartition & domain ) const;

  void computeTolerances( DomainPartition & domain ) const;

  void computeFaceNodalArea( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                             ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                             localIndex const kf0,
                             array1d< real64 > & nodalArea ) const;

  void computeProjectors( localIndex const faceIndex,
                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoords,
                          ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                          ArrayOfArraysView< localIndex const > const & faceToEdgeMap,
                          arrayView2d< localIndex const > const & edgeToNodeMap,
                          arrayView2d< real64 const > const faceCenters,
                          arrayView2d< real64 const > const faceNormals,
                          arrayView1d< real64 const > const faceAreas,
                          array1d< real64 > & basisIntegrals );                           

  void computeFaceIntegrals( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoords,
                      localIndex const (&faceToNodes)[11],
                      localIndex const (&faceToEdges)[11],
                      localIndex const & numFaceVertices,
                      real64 const & faceArea,
                      real64 const (&faceCenter)[3],
                      real64 const (&faceNormal)[3],
                      arrayView2d< localIndex const > const & edgeToNodes,
                      real64 const & invCellDiameter,
                      real64 const (&cellCenter)[3],
                      array1d< real64 > & basisIntegrals,
                      real64 (& threeDMonomialIntegrals)[3] );
  
  real64 const machinePrecision = std::numeric_limits< real64 >::epsilon();

  string getStabilizationName() const { return m_stabilizationName; }

protected:
  virtual void postProcessInput() override final;

private:
  string m_stabilizationName;

  string m_polyMeshInfSupName; // Testing the flexibility of the polyhedral mesh

  localIndex m_contactRelationFullIndex;

  real64 const m_slidingCheckTolerance = 0.05;

  real64 m_initialResidual[3] = {0.0, 0.0, 0.0};

  static const localIndex m_MFN; // Maximum number of nodes on a contact face

  void createPreconditioner( DomainPartition const & domain );

  void computeFaceDisplacementJump( DomainPartition & domain );

  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;

  template< localIndex DIMENSION, typename POINT_COORDS_TYPE >
  GEOS_HOST_DEVICE
  inline static real64 computeDiameter( POINT_COORDS_TYPE points,
                                localIndex const & numPoints )
  {
    real64 diameter = 0;
    for( localIndex numPoint = 0; numPoint < numPoints; ++numPoint )
    {
      for( localIndex numOthPoint = 0; numOthPoint < numPoint; ++numOthPoint )
      {
        real64 candidateDiameter = 0.0;
        for( localIndex i = 0; i < DIMENSION; ++i )
        {
          real64 coordDiff = points[numPoint][i] - points[numOthPoint][i];
          candidateDiameter += coordDiff * coordDiff;
        }
        if( diameter < candidateDiameter )
        {
          diameter = candidateDiameter;
        }
      }
    }
    return LvArray::math::sqrt< real64 >( diameter );
  }

  struct viewKeyStruct : ContactSolverBase::viewKeyStruct
  {
    constexpr static char const * stabilizationNameString() { return "stabilizationName"; }
    constexpr static char const * contactRelationNameString() { return "contactRelationName"; }
    constexpr static char const * activeSetMaxIterString() { return "activeSetMaxIter"; } // TODO: remove
    
    constexpr static char const * rotationMatrixString() { return "rotationMatrix"; }

    constexpr static char const * slidingCheckToleranceString() { return "slidingCheckTolerance"; }
    constexpr static char const * normalDisplacementToleranceString() { return "normalDisplacementTolerance"; }
    constexpr static char const * normalTractionToleranceString() { return "normalTractionTolerance"; }
    constexpr static char const * slidingToleranceString() { return "slidingTolerance"; }
    constexpr static char const * usePolyMeshForInfSupName() { return "polyMeshInfSup"; } // Testing the flexibility of the polyhedral mesh

    static constexpr char const * transMultiplierString() { return "penaltyStiffnessTransMultiplier"; }

  };

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSCONFORMINGFRACTURESVEM_HPP_ */
