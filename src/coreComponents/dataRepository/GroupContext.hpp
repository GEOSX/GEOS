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
 * @file GroupContext.hpp
 */

#ifndef GEOS_DATAREPOSITORY_GROUPCONTEXT_HPP_
#define GEOS_DATAREPOSITORY_GROUPCONTEXT_HPP_

#include "DataContext.hpp"
#include "Group.hpp"

namespace geos
{
namespace dataRepository
{


/**
 * @class GroupContext
 *
 * Helps to know where a Group is in the hierarchy.
 * See DataContext class for more info.
 */
class GroupContext : public DataContext
{
public:

  /**
   * @brief Construct a new GroupContext object
   * @param group The reference to the Group related to this GroupContext.
   */
  GroupContext( Group & group );

  /**
   * @return the reference to the Group related to this GroupContext.
   */
  Group const & getGroup() const;

protected:

  /**
   * @brief Construct a new GroupContext object
   * @param group The reference to the Group related to this GroupContext.
   * @param objectName Target object name.
   */
  GroupContext( Group & group, string const & objectName );

  /// The reference to the Group related to this GroupContext.
  Group & m_group;

private:

  /**
   * @return the group path with the file & line of the first parent for which this information exists.
   */
  string toString() const override;
  /**
   * @copydoc DataContext::getToStringInfo()
   */
  ToStringInfo getToStringInfo() const override;
};


} /* namespace dataRepository */
} /* namespace geos */

#endif /* GEOS_DATAREPOSITORY_GROUPCONTEXT_HPP_ */
