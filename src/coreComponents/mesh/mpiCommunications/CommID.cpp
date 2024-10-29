/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "CommID.hpp"

#include "common/Logger.hpp"

#include <set>

namespace geos
{

CommID::CommID( std::set< int > & freeIDs ):
  m_freeIDs( freeIDs ),
  m_id( -1 )
{
  GEOS_ERROR_IF_EQ( freeIDs.size(), 0 );
  m_id = *freeIDs.begin();
  freeIDs.erase( freeIDs.begin() );
}

CommID::CommID( CommID && src ):
  m_freeIDs( src.m_freeIDs ),
  m_id( src.m_id )
{
  src.m_id = -1;
}

CommID::~CommID()
{
  GEOS_ERROR_IF( m_freeIDs.count( m_id ) > 0, "Attempting to release commID that is already free: " << m_id );

  m_freeIDs.insert( m_id );
  m_id = -1;
}

} /* namespace geos */
