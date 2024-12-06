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
 * @file CatalystActions.hpp
 */

// source includes
#include "ConduitCapsule.hpp"

namespace geos
{

/**
 * @brief Wrapped call to catalyst_initialize
 * @param the encapsulated conduit node to convert and send to initialize method
 * @returns true on success
 */
bool CatalystInitialize(ConduitCapsule* simulationNode);

/**
 * @brief Wrapped call to catalyst_execute
 * @param the encapsulated conduit node to convert and send to execute method
 * @returns true on success
 */
bool CatalystExecute(ConduitCapsule* simulationNode);

/**
 * @brief Wrapped call to catalyst_finalize
 * @param the encapsulated conduit node to convert and send to finalize method
 * @returns true on success
 */
bool CatalystFinalize(ConduitCapsule* simulationNode);

}
