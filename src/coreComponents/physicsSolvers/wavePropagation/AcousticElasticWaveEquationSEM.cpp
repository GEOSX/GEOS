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
 * @file AcousticElasticWaveEquation.cpp
 */

#include "AcousticElasticWaveEquationSEM.hpp"
#include "dataRepository/Group.hpp"

namespace geos
{
using namespace dataRepository;
using namespace fields;

AcousticElasticWaveEquationSEM::AcousticElasticWaveEquationSEM( const string & name,
                                                                Group * const parent )
  : Base( name, parent )
{
  // std::cout << "\t[AcousticElasticWaveEquationSEM::AcousticElasticWaveEquationSEM]" << std::endl;
}

void AcousticElasticWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups();

  auto acousNodesSet = acousticSolver()->getSolverNodesSet();
  auto elasNodesSet = elasticSolver()->getSolverNodesSet();

  for (auto val : acousNodesSet)
  {
    if (elasNodesSet.contains(val)) m_interfaceNodesSet.insert(val);
  }

  // std::cout << "\t[AcousticElasticWaveEquationSEM::initializePostInitialConditionsPreSubGroups] size=" << m_solverTargetNodesSet.size() << std::endl;
  GEOS_THROW_IF( m_interfaceNodesSet.size() == 0, "Failed to compute interface: check xml input (solver order)", std::runtime_error );
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticElasticWaveEquationSEM, string const &, Group * const )

} /* namespace geos */
