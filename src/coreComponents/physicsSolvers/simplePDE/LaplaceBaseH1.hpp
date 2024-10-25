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

#ifndef GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_BASE_HPP
#define GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_BASE_HPP

#include "codingUtilities/EnumStrings.hpp"   // facilities for enum-string conversion (for reading enum values from XML input)
#include "physicsSolvers/SolverBase.hpp"  // an abstraction class shared by all physics solvers
#include "fieldSpecification/FieldSpecificationManager.hpp" // a manager that can access and set values on the discretized domain

namespace geos
{

// Like most physics solvers, the Laplace solver derives from a generic SolverBase class.
// The base class is densely Doxygen-commented and worth a look if you have not done so already.
// Most important system assembly steps, linear and non-linear resolutions, and time-stepping mechanisms
// are implemented at the SolverBase class level and can thus be used in Laplace without needing reimplementation.
class LaplaceBaseH1 : public SolverBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  LaplaceBaseH1() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  LaplaceBaseH1( const string & name,
                 Group * const parent );

  /// Destructor
  virtual ~LaplaceBaseH1() override;

//START_SPHINX_INCLUDE_REGISTERDATAONMESH

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override final;

//END_SPHINX_INCLUDE_REGISTERDATAONMESH

  //START_SPHINX_INCLUDE_SOLVERINTERFACE
  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  applyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override final;

  virtual void
    resetStateToBeginningOfStep( DomainPartition & GEOS_UNUSED_PARAM( domain ) ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  /// This method is specific to this Laplace solver.
  /// It is used to apply Dirichlet boundary condition
  /// and called when the base class applyBoundaryConditions() is called.
  virtual void applyDirichletBCImplicit( real64 const time,
                                         DofManager const & dofManager,
                                         DomainPartition & domain,
                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                         arrayView1d< real64 > const & localRhs );

  //END_SPHINX_INCLUDE_SOLVERINTERFACE

  /// Choice of transient treatment options (steady or backward Euler scheme):
  //START_SPHINX_INCLUDE_TIMEINTOPT
  enum class TimeIntegrationOption : integer
  {
    SteadyState,
    ImplicitTransient
  };
  //END_SPHINX_INCLUDE_TIMEINTOPT


  /// This structure stores ``dataRepository::ViewKey`` objects used as binding between the input
  /// XML tags and source code variables (here, timeIntegrationOption and fieldVarName)
  //START_SPHINX_INCLUDE_VIEWKEY
  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    static constexpr char const * timeIntegrationOption() { return "timeIntegrationOption"; }
    static constexpr char const * fieldVarName() { return "fieldName"; }
  };
  //END_SPHINX_INCLUDE_VIEWKEY

protected:

  // These two data members are specific to the Laplace solver:
  /// User-defined name of the physical quantity we wish to solve for (such as "Temperature", etc.)
  string m_fieldName;
  /// Choice of transient treatment (SteadyState, ImplicitTransient)
  TimeIntegrationOption m_timeIntegrationOption;

};


/* REGISTERING NEW ENUMS:
   --------------------------
   In order to register an enumeration type with the Data Repository and have its value read from input,
   we must define stream insertion/extraction operators. This is a common task, so GEOSX provides
   a facility for automating it. Upon including ``common/EnumStrings.hpp``, we can call the following macro
   at the namespace scope (in this case, right after the ``LaplaceBaseH1`` class definition is complete):
 */
//START_SPHINX_INCLUDE_REGENUM
ENUM_STRINGS( LaplaceBaseH1::TimeIntegrationOption,
              "SteadyState",
              "ImplicitTransient" );
//END_SPHINX_INCLUDE_REGENUM

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_BASE_HPP */
