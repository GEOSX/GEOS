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
 * @file HydrofractureSolver.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVER_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "dataRepository/LogLevelsInfo.hpp"

namespace geos
{

using dataRepository::Group;

template< typename POROMECHANICS_SOLVER = SinglePhasePoromechanics<> >
class HydrofractureSolver : public POROMECHANICS_SOLVER
{
public:

  using Base = POROMECHANICS_SOLVER;
  using Base::m_solvers;
  using Base::m_names;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;
  using Base::m_linearSolverParameters;

  using Base::registerWrapper;
  using Base::forDiscretizationOnMeshTargets;
  using Base::getMeshModificationTimestamp;
  using Base::getSystemSetupTimestamp;
  using Base::nonlinearImplicitStep;
  using Base::implicitStepComplete;
  using Base::getLogLevel;
  using Base::setSystemSetupTimestamp;
  using Base::setupDofs;
  using Base::flowSolver;
  using Base::solidMechanicsSolver;
  using Base::assembleElementBasedTerms;


  /**
   * @brief main constructor for HydrofractureSolver objects
   * @param name the name of this instantiation of HydrofractureSolver in the repository
   * @param parent the parent group of this instantiation of HydrofractureSolver
   */
  HydrofractureSolver( const string & name,
                       Group * const parent );

  /// Destructor for the class
  ~HydrofractureSolver() override {}

  static string catalogName()
  {
    // single phase
    if constexpr ( std::is_same_v< POROMECHANICS_SOLVER, SinglePhasePoromechanics< SinglePhaseBase > > )
    {
      return "Hydrofracture";
    }
//  // multi phase (TODO)
//  else if constexpr ( std::is_same_v< POROMECHANICS_SOLVER, MultiphasePoromechanics< CompositionalMultiphaseBase > > )
//  {
//    return "MultiphaseHydrofracture";
//  }
  }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "poromechanics"; }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void registerDataOnMesh( Group & MeshBodies ) override final;

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override final;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain ) override final;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual real64 setNextDt( real64 const & currentTime,
                            real64 const & lastDt,
                            DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override final;

  virtual void implicitStepComplete( real64 const & time_n,
                                     real64 const & dt,
                                     DomainPartition & domain ) override final;

  /**@}*/

  void updateHydraulicApertureAndFracturePermeability( DomainPartition & domain );

  void assembleForceResidualDerivativeWrtPressure( DomainPartition & domain,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs );

  void assembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const & domain,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix );

  std::unique_ptr< CRSMatrix< real64, localIndex > > & getRefDerivativeFluxResidual_dAperture()
  {
    return m_derivativeFluxResidual_dAperture;
  }

  CRSMatrixView< real64, localIndex const > getDerivativeFluxResidual_dNormalJump()
  {
    return m_derivativeFluxResidual_dAperture->toViewConstSizes();
  }

  CRSMatrixView< real64 const, localIndex const > getDerivativeFluxResidual_dNormalJump() const
  {
    return m_derivativeFluxResidual_dAperture->toViewConst();
  }

  enum class InitializationType : integer
  {
    Pressure,
    Displacement,
  };

  struct viewKeyStruct : Base::viewKeyStruct
  {
    constexpr static char const * surfaceGeneratorNameString() { return "surfaceGeneratorName"; }

    constexpr static char const * maxNumResolvesString() { return "maxNumResolves"; }

    constexpr static char const * isMatrixPoroelasticString() { return "isMatrixPoroelastic"; }

    constexpr static char const * newFractureInitializationTypeString() { return "newFractureInitializationType"; }

    constexpr static char const * useQuasiNewtonString() { return "useQuasiNewton"; }

    constexpr static char const * isLaggingFractureStencilWeightsUpdateString() { return "isLaggingFractureStencilWeightsUpdate"; }

#ifdef GEOS_USE_SEPARATION_COEFFICIENT
    constexpr static char const * separationCoeff0String() { return "separationCoeff0"; }
    constexpr static char const * apertureAtFailureString() { return "apertureAtFailure"; }
#endif
  };

protected:

  virtual void postInputInitialization() override final;

  /**
   * @Brief add the nnz induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param rowLenghts the nnz in each row
   */
  void addFluxApertureCouplingNNZ( DomainPartition & domain,
                                   DofManager & dofManager,
                                   arrayView1d< localIndex > const & rowLengths ) const;


  /**
   * @Brief add the sparsity pattern induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  void addFluxApertureCouplingSparsityPattern( DomainPartition & domain,
                                               DofManager & dofManager,
                                               SparsityPatternView< globalIndex > const & pattern ) const;


  void setUpDflux_dApertureMatrix( DomainPartition & domain,
                                   DofManager const & dofManager,
                                   CRSMatrix< real64, globalIndex > & localMatrix );


private:

  virtual real64 fullyCoupledSolverStep( real64 const & time_n,
                                         real64 const & dt,
                                         int const cycleNumber,
                                         DomainPartition & domain ) override final;


  /**
   * @brief Initialize fields on the newly created elements of the fracture.
   * @param domain the physical domain object
   */
  void initializeNewFractureFields( DomainPartition & domain );

  // name of the contact relation
  string m_contactRelationName;

  /// name of the surface generator
  string m_surfaceGeneratorName;

  /// pointer to the surface generator
  SurfaceGenerator * m_surfaceGenerator;

  // it is only important for this case.
  std::unique_ptr< CRSMatrix< real64, localIndex > > m_derivativeFluxResidual_dAperture;

  integer m_maxNumResolves;
  integer m_numResolves[2];

  integer m_isMatrixPoroelastic;

  // flag to determine which initialization type to use for the new fracture cell
  InitializationType m_newFractureInitializationType;

  integer m_useQuasiNewton;   // use Quasi-Newton (see https://arxiv.org/abs/2111.00264)

  // flag to determine whether or not to apply lagging update for the fracture stencil weights
  integer m_isLaggingFractureStencilWeightsUpdate;

};

ENUM_STRINGS( HydrofractureSolver< SinglePhasePoromechanics< SinglePhaseBase > >::InitializationType,
              "Pressure",
              "Displacement" );


} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVER_HPP_ */
