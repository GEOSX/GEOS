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
 * @file SinglePhaseReservoirAndWells.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASERESERVOIRANDWELLS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASERESERVOIRANDWELLS_HPP_

#include "physicsSolvers/multiphysics/CoupledReservoirAndWellsBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"

namespace geos
{

/// @tparam RESERVOIR_SOLVER single-phase flow or single-phase poromechanics solver
template< typename RESERVOIR_SOLVER = SinglePhaseBase >
class SinglePhaseReservoirAndWells : public CoupledReservoirAndWellsBase< RESERVOIR_SOLVER,
                                                                          SinglePhaseWell >
{
public:

  using Base = CoupledReservoirAndWellsBase< RESERVOIR_SOLVER,
                                             SinglePhaseWell >;
  using Base::m_solvers;
  using Base::m_linearSolverParameters;

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  SinglePhaseReservoirAndWells( const string & name,
                                dataRepository::Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseReservoirAndWells() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    if constexpr (std::is_same_v< RESERVOIR_SOLVER, SinglePhaseBase > ) // special case
    {
      return "SinglePhaseReservoir";
    }
    else // default
    {
      return RESERVOIR_SOLVER::catalogName() + "Reservoir";
    }
  }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void addCouplingSparsityPattern( DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           SparsityPatternView< globalIndex > const & pattern ) const override;

  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override;

  void
  assembleHydrofracFluxTerms( real64 const time_n,
                              real64 const dt,
                              DomainPartition const & domain,
                              DofManager const & dofManager,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs,
                              CRSMatrixView< real64, localIndex const > const & dR_dAper )
  { flowSolver()->assembleHydrofracFluxTerms( time_n, dt, domain, dofManager, localMatrix, localRhs, dR_dAper ); }

  template< typename SUBREGION_TYPE >
  void accumulationAssemblyLaunch( DofManager const & dofManager,
                                   SUBREGION_TYPE const & subRegion,
                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                   arrayView1d< real64 > const & localRhs )
  { flowSolver()->accumulationAssemblyLaunch( dofManager, subRegion, localMatrix, localRhs ); }

  void prepareStencilWeights( DomainPartition & domain ) const
  { flowSolver()->prepareStencilWeights( domain ); }
  void updateStencilWeights( DomainPartition & domain ) const
  { flowSolver()->updateStencilWeights( domain ); }

protected:

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

private:

  SinglePhaseBase * flowSolver() const;

  void setMGRStrategy();

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASERESERVOIRANDWELLS_HPP_ */
