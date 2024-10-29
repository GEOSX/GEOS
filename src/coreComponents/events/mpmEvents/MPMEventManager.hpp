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


#ifndef GEOS_MPMEVENTMANAGER_HPP_
#define GEOS_MPMEVENTMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "MPMEventBase.hpp"

namespace geos
{

namespace dataRepository
{
namespace keys
{
string const MPMEvents( "MPMEvents" );
}
}

/**
 * @class MPMEventManager
 *
 * A class for managing mpm events.
 */
class MPMEventManager : public dataRepository::Group
{
public:
  /**
   * @brief Constructor for the MPMEventManager
   * @param[in] name the name of the MPMEventManager
   * @param[in] parent group this MPMEventManager
   */
  MPMEventManager( string const & name,
                Group * const parent );

  /**
   * @brief Default destructor for the EventManager
   */
  virtual ~MPMEventManager() override;

  /**
   * @brief Create a child MPMEventManager
   * @param[in] childKey the key of the Event to be added
   * @param[in] childName the name of the child to be added
   * @return the Event
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief This method is used to expand any catalogs in the data structure
   */
  virtual void expandObjectCatalogs() override;

  /**
   * @name viewKeyStruct/groupKeyStruct
   */
  ///@{
  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    // static constexpr char const * minTimeString() { return "minTime"; }
    // dataRepository::ViewKey minTime = { "minTime" };
  } viewKeys;
  /// @endcond
  ///@}

  /// Alias to access the object catalog for EventBase derived types.
  using CatalogInterface = dataRepository::CatalogInterface< MPMEventBase, string const &, Group * const >;

  /// @copydoc dataRepository::Group::getCatalog()
  static CatalogInterface::CatalogType & getCatalog();
};


} /* namespace geos */

#endif /* GEOS_MPMEVENTMANAGER_HPP_ */
