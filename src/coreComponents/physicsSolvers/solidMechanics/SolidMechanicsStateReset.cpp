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
 * @file SolidMechanicsStateReset.cpp
 */

#include "SolidMechanicsStateReset.hpp"

#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

SolidMechanicsStateReset::SolidMechanicsStateReset( const string & name,
                                                    Group * const parent ):
  TaskBase( name, parent ),
  m_solidSolverName()
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver" );

  registerWrapper( viewKeyStruct::resetDisplacementsString(), &m_resetDisplacements ).
    setApplyDefaultValue( true ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to reset displacements (and velocities)" );

  registerWrapper( viewKeyStruct::disableInelasticityString(), &m_disableInelasticity ).
    setApplyDefaultValue( false ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to enable/disable inelastic behavior" );

  addLogLevel< logInfo::Initialization >();
}

SolidMechanicsStateReset::~SolidMechanicsStateReset()
{}

void SolidMechanicsStateReset::postInputInitialization()
{
  Group & problemManager = this->getGroupByPath( "/Problem" );
  Group & physicsSolverManager = problemManager.getGroup( "Solvers" );

  GEOS_THROW_IF( !physicsSolverManager.hasGroup( m_solidSolverName ),
                 GEOS_FMT( "Task {}: physics solver named {} not found",
                           getDataContext(), m_solidSolverName ),
                 InputError );

  m_solidSolver = &physicsSolverManager.getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
}

bool SolidMechanicsStateReset::execute( real64 const time_n,
                                        real64 const GEOS_UNUSED_PARAM( dt ),
                                        integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                        integer const GEOS_UNUSED_PARAM( eventCounter ),
                                        real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                        DomainPartition & domain )
{
  m_solidSolver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                               MeshLevel & mesh,
                                                                               arrayView1d< string const > const & regionNames )
  {
    // Option 1: zero out velocity, incremental displacement, and displacement
    if( m_resetDisplacements )
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Initialization,
                                  GEOS_FMT( "Task `{}`: at time {}s, physics solver `{}` is resetting total displacement and velocity to zero",
                                            getName(), time_n, m_solidSolverName ) );

      NodeManager & nodeManager = mesh.getNodeManager();

      if( nodeManager.hasField< solidMechanics::velocity >() )
      {
        nodeManager.getField< solidMechanics::velocity >().zero();
      }
      nodeManager.getField< solidMechanics::totalDisplacement >().zero();
      nodeManager.getField< solidMechanics::incrementalDisplacement >().zero();
    }

    // Option 2: enable / disable inelastic behavior
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            CellElementSubRegion & subRegion )
    {
      string const & solidMaterialName = subRegion.getReference< string >( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString() );
      Group & constitutiveModels = subRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );

      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Initialization, GEOS_FMT( "Task `{}`: at time {}s, solid model `{}` is setting inelastic behavior to `{}` on subRegion `{}`. ",
                                                                     getName(), time_n, solidMaterialName,
                                                                     m_disableInelasticity ? "OFF" : "ON",
                                                                     subRegion.getName() ) );

      SolidBase & constitutiveRelation = constitutiveModels.getGroup< SolidBase >( solidMaterialName );
      constitutiveRelation.disableInelasticity( m_disableInelasticity );
    } );
  } );

  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        SolidMechanicsStateReset,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
