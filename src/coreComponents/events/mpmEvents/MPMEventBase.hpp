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
 * @file MPMEventBase.hpp
 */

#ifndef GEOSX_MPMEVENTBASE_HPP_
#define GEOSX_MPMEVENTBASE_HPP_

#include "dataRepository/Group.hpp"

namespace geos
{

/**
 * @class MPMEventBase
 *
 * This class implements the abstract base mpm event for the solid mechanics material point method solver
 */
class MPMEventBase : public dataRepository::Group
{
public:
  /**
   * Constructor
   * @param time at which event starts
   * @param interval time interval between during which the event is performed
   */
  MPMEventBase( string const & name,
                Group * const parent );

  /// Destructor
  virtual ~MPMEventBase() override;


  /**
   * @copydoc dataRepository::Group::createChild()
   *
   * An event may have an arbitrary number of sub-events defined as children in the input xml.
   * e.g.:
   * @code{.unparsed}
   * <Events>
   *   <PeriodicEvent name="base_event" ...>
   *     <PeriodicEvent name="sub_event" .../>
   *     ...
   *   </PeriodicEvent>
   * </Events>
   * @endcode
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;


  /**
   * @brief Expand any catalogs in the data structure.
   */
  virtual void expandObjectCatalogs() override;


  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "MPMEventBase"; }

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * timeString() { return "time"; }
    static constexpr char const * intervalString() { return "interval"; }
    static constexpr char const * isCompleteString() { return "isComplete"; }

    dataRepository::ViewKey time = { timeString() };
    dataRepository::ViewKey interval = { intervalString() };
    dataRepository::ViewKey isComplete = { isCompleteString() };
  } viewKeys;
  /// @endcond


  /// Catalog interface
  using CatalogInterface = dataRepository::CatalogInterface< MPMEventBase, string const &, Group * const >;
  /// @copydoc dataRepository::Group::getCatalog()
  static CatalogInterface::CatalogType & getCatalog();

  real64 getTime() const { return m_time; }
  real64 getInterval() const { return m_interval; }

  int isComplete() const { return m_isComplete; }
  void setIsComplete( int isComplete ) { m_isComplete = isComplete; }

protected:
  // Event variables
  real64 m_time;
  real64 m_interval;
  int m_isComplete;
};

} /* namespace geos */

#endif /* GEOSX_MPMEVENTBASE_HPP_ */
