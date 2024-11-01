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
 * @file OutputManager.hpp
 */

#ifndef GEOS_FILEIO_OUTPUTS_OUTPUTMANAGER_HPP_
#define GEOS_FILEIO_OUTPUTS_OUTPUTMANAGER_HPP_

#include "dataRepository/Group.hpp"


namespace geos
{
namespace dataRepository
{
namespace keys
{}
}

/**
 * @class OutputManager
 *
 * An class for managing output types
 */
class OutputManager : public dataRepository::Group
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  OutputManager( string const & name,
                 Group * const parent );

  /// Destructor
  virtual ~OutputManager() override;

  /// @copydoc geos::dataRepository::Group::createChild( string const & childKey, string const & childName )
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void expandObjectCatalogs() override;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    dataRepository::ViewKey time = { "time" };
  } viewKeys;
  /// @endcond
};


} /* namespace geos */

#endif /* GEOS_FILEIO_OUTPUTS_OUTPUTMANAGER_HPP_ */
