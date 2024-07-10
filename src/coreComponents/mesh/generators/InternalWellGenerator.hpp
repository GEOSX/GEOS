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

/*
 * @file InternalWellGenerator.hpp
 *
 */

#ifndef GEOS_MESH_GENERATORS_INTERNALWELLGENERATOR_HPP_
#define GEOS_MESH_GENERATORS_INTERNALWELLGENERATOR_HPP_

#include "WellGeneratorBase.hpp"

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"


namespace geos
{

/**
 * @class InternalWellGenerator
 *
 * This class processes the data of a single well from the XML and generates the well geometry
 */
class InternalWellGenerator : public WellGeneratorBase
{
public:


  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  InternalWellGenerator( const string & name,
                         Group * const parent );

  /**
   * @brief Get the catalog name.
   * @return the name of this type in the catalog
   */
  static string catalogName() { return "InternalWell"; }

  ///@}

protected:
  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void postInputInitialization() override;

};
}
#endif /* GEOS_MESH_GENERATORS_INTERNALWELLGENERATOR_HPP_ */
