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

#include "mesh/generators/LineBlock.hpp"


namespace geos
{

LineBlock::LineBlock( string const & name, Group * const parent ):
  LineBlockABC( name, parent )
{}

void LineBlock::setPrevElemIndices( arrayView1d< arrayView1d< globalIndex const > const > prevElemIndices )
{
  int size = prevElemIndices.size();
  m_prevElemId.resize( size );
  for( int i = 0; i < size; i++ )
  {
    m_prevElemId[i] = prevElemIndices[i];
  }
}
}
