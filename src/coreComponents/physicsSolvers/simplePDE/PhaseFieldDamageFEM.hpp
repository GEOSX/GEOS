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
 * @file PhaseFieldDamageFEM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGE_HPP_
#define GEOS_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGE_HPP_

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/SolverBase.hpp"

struct stabledt
{
  double m_maxdt;
};

namespace geos
{
namespace dataRepository
{
class Group;
}
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;

class PhaseFieldDamageFEM : public SolverBase
{
public:
  PhaseFieldDamageFEM( const string & name, Group * const parent );

  virtual ~PhaseFieldDamageFEM() override;

  static string catalogName()
  {
    return "PhaseFieldDamageFEM";
  }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  static string coupledSolverAttributePrefix() { return "damage"; }

  virtual void registerDataOnMesh( Group & meshBodies ) override final;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived
   * classes
   */
  /**@{*/

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  virtual real64 explicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain ) override;

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void assembleSystem( real64 const time, real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void applyBoundaryConditions( real64 const time, real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override;

  virtual real64 calculateResidualNorm( real64 const & time_n,
                                        real64 const & dt,
                                        DomainPartition const & domain,
                                        DofManager const & dofManager,
                                        arrayView1d< real64 const > const & localRhs ) override;

  virtual void applySystemSolution( DofManager const & dofManager,
                                    arrayView1d< real64 const > const & localSolution,
                                    real64 const scalingFactor,
                                    real64 const dt,
                                    DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override final;

  virtual void
  implicitStepSetup( real64 const &,
                     real64 const &,
                     DomainPartition & ) override {}

  virtual void implicitStepComplete( real64 const & time,
                                     real64 const & dt,
                                     DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & ) override {}

  /**@}*/

  void applyDirichletBCImplicit( real64 const time,
                                 DofManager const & dofManager,
                                 DomainPartition & domain,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs );

  void applyIrreversibilityConstraint( DofManager const & dofManager,
                                       DomainPartition & domain,
                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                       arrayView1d< real64 > const & localRhs );

  virtual void saveSequentialIterationState( DomainPartition & domain ) override;

  enum class TimeIntegrationOption
  {
    SteadyState,
    ImplicitTransient,
    ExplicitTransient
  };

  enum class LocalDissipation
  {
    Linear,
    Quadratic,
  };

  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    static constexpr char const * coeffNameString() { return "coeffField"; }
    static constexpr char const * localDissipationOptionString() { return "localDissipation"; }
    static constexpr char const * irreversibilityFlagString() { return "irreversibilityFlag"; }
    static constexpr char const * damageUpperBoundString() { return "damageUpperBound"; }
    static constexpr char const * fracturePressureTermFlagString() { return "fracturePressureTermFlag"; }
    static constexpr char const * solidModelNamesString() { return "solidMaterialNames"; }

    dataRepository::ViewKey timeIntegrationOption = { "timeIntegrationOption" };
    dataRepository::ViewKey fieldVarName = { "fieldName" };
  } PhaseFieldDamageFEMViewKeys;

  inline ParallelVector const * getSolution() const
  {
    return &m_solution;
  }

  inline globalIndex getSize() const
  {
    return m_matrix.numGlobalRows();
  }

  string const & getFieldName() const
  {
    return m_fieldName;
  }

protected:
  virtual void postInputInitialization() override final;

private:
  string m_fieldName;
  stabledt m_stabledt;
  TimeIntegrationOption m_timeIntegrationOption;
  LocalDissipation m_localDissipationOption;
  integer m_irreversibilityFlag;
  real64 m_damageUpperBound;
  integer m_fracturePressureTermFlag;

  array1d< real64 > m_coeff;

  PhaseFieldDamageFEM();
};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( PhaseFieldDamageFEM::LocalDissipation,
              "Linear",
              "Quadratic" );
ENUM_STRINGS( PhaseFieldDamageFEM::TimeIntegrationOption,
              "SteadyState",
              "ImplicitTransient",
              "ExplicitTransient" );

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGE_HPP_ */
