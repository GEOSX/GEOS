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
 * @file PhaseFieldFractureSolver.cpp
 *
 */

#include "PhaseFieldFractureSolver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

PhaseFieldFractureSolver::PhaseFieldFractureSolver( const string & name,
                                                    Group * const parent ):
  Base( name, parent )
{}

PhaseFieldFractureSolver::~PhaseFieldFractureSolver()
{
  // TODO Auto-generated destructor stub
}

void PhaseFieldFractureSolver::postInputInitialization()
{
  Base::postInputInitialization();
  GEOS_WARNING_IF( getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::FullyImplicit,
                   "FullyImplicit coupling not implemented for this solver. A sequential coupling approach will be used." );
  getNonlinearSolverParameters().m_couplingType = NonlinearSolverParameters::CouplingType::Sequential;
}

void PhaseFieldFractureSolver::mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType )
{

  GEOS_MARK_FUNCTION;
  if( solverType ==  static_cast< integer >( SolverType::Damage ) )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & regionNames )
    {
      NodeManager & nodeManager = mesh.getNodeManager();

      string const & damageFieldName = damageSolver()->getFieldName();

      string const & discretizationName = damageSolver()->getDiscretizationName();

      //should get reference to damage field here.
      arrayView1d< real64 const > const nodalDamage = nodeManager.getReference< array1d< real64 > >( damageFieldName );

      ElementRegionManager & elemManager = mesh.getElemManager();

      // begin region loop
      elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [discretizationName, nodalDamage]
                                                                  ( localIndex const,
                                                                  CellElementSubRegion & elementSubRegion )
      {
        string const & solidModelName = elementSubRegion.getReference< string >( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString());
        constitutive::SolidBase &
        solidModel = elementSubRegion.getConstitutiveModel< constitutive::SolidBase >( solidModelName );

        ConstitutivePassThru< DamageBase >::execute( solidModel, [&elementSubRegion, discretizationName, nodalDamage]( auto & damageModel )
        {
          using CONSTITUTIVE_TYPE = TYPEOFREF( damageModel );
          typename CONSTITUTIVE_TYPE::KernelWrapper constitutiveUpdate = damageModel.createKernelUpdates();

          arrayView2d< real64 > const damageFieldOnMaterial = constitutiveUpdate.m_damage;
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodes = elementSubRegion.nodeList();

          finiteElement::FiniteElementBase const &
          fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( discretizationName );

          finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( fe, [=, &elementSubRegion] ( auto & finiteElement )
          {
            using FE_TYPE = TYPEOFREF( finiteElement );

            DamageInterpolationKernel< FE_TYPE > interpolationKernel( elementSubRegion );

            interpolationKernel.interpolateDamage( elemToNodes, nodalDamage, damageFieldOnMaterial );
          } );
        } );
      } );
    } );
  }
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, PhaseFieldFractureSolver, string const &, Group * const )

} /* namespace geos */
