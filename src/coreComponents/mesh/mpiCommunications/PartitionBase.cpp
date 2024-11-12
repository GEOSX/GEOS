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

#include "PartitionBase.hpp"

namespace geos
{

PartitionBase::PartitionBase()
  :
  m_rank( MpiWrapper::commRank( MPI_COMM_GEOS ) ),
  m_size( MpiWrapper::commSize( MPI_COMM_GEOS ) ),
  m_numColors( -1 ),
  m_color( -1 ),
  m_min{ 0.0, 0.0, 0.0 },
  m_max{ 0.0, 0.0, 0.0 },
  m_gridSize{ 0.0, 0.0, 0.0 },
  m_gridMin{ 0.0, 0.0, 0.0 },
  m_gridMax{ 0.0, 0.0, 0.0 }
{}

PartitionBase::~PartitionBase()
{}

}
