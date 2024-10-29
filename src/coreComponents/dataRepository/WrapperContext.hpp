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

/**
 * @file WrapperContext.hpp
 */

#ifndef GEOS_DATAREPOSITORY_WRAPPERCONTEXT_HPP_
#define GEOS_DATAREPOSITORY_WRAPPERCONTEXT_HPP_

#include "GroupContext.hpp"
#include "WrapperBase.hpp"

namespace geos
{
namespace dataRepository
{


/**
 * @class WrapperContext
 *
 * Dedicated implementation of GroupContext for Wrapper.
 * See DataContext class for more info.
 */
class WrapperContext final : public GroupContext
{
public:

  /**
   * @brief Construct a new WrapperContext object
   * @param wrapper the target Wrapper object
   */
  WrapperContext( WrapperBase & wrapper );

private:

  string const m_typeName;

  /**
   * @return the parent group DataContext followed by the wrapper name.
   */
  string toString() const override;

};


} /* namespace dataRepository */
} /* namespace geos */

#endif /* GEOS_DATAREPOSITORY_WRAPPERCONTEXT_HPP_ */
