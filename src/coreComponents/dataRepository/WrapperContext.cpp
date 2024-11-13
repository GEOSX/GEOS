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
 * @file WrapperContext.cpp
 */

#include "WrapperContext.hpp"

namespace geos
{
namespace dataRepository
{


WrapperContext::WrapperContext( WrapperBase & wrapper ):
  GroupContext( wrapper.getParent(), wrapper.getParent().getName() + '/' + wrapper.getName() ),
  m_typeName( wrapper.getName() )
{}

string WrapperContext::toString() const
{
  ToStringInfo const info = m_group.getDataContext().getToStringInfo();
  return info.hasInputFileInfo() ?
         GEOS_FMT( "{} ({}, l.{})", m_targetName, info.m_filePath, info.m_line ) :
         GEOS_FMT( "{}/{}", m_group.getDataContext().toString(), m_typeName );
}


} /* namespace dataRepository */
} /* namespace geos */
