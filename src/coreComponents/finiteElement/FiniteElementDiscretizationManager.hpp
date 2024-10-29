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
 * @file FiniteElementSpaceManager.hpp
 */

#ifndef GEOS_FINITEELEMENT_FINITEELEMENTSPACEMANAGER_HPP_
#define GEOS_FINITEELEMENT_FINITEELEMENTSPACEMANAGER_HPP_

#include "dataRepository/Group.hpp"

namespace geos
{

class FiniteElementDiscretizationManager : public dataRepository::Group
{
public:
  FiniteElementDiscretizationManager() = delete;
  FiniteElementDiscretizationManager( string const & name, Group * const parent );
  virtual ~FiniteElementDiscretizationManager() override;

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void expandObjectCatalogs() override;

};

} /* namespace geos */

#endif /* GEOS_FINITEELEMENT_FINITEELEMENTSPACEMANAGER_HPP_ */
