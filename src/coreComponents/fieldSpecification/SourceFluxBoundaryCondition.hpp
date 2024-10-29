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

/*
 * SourceFluxBoundaryCondition.hpp
 *
 */

#ifndef GEOS_FIELDSPECIFICATION_SOURCEFLUXBOUNDARYCONDITION_HPP_
#define GEOS_FIELDSPECIFICATION_SOURCEFLUXBOUNDARYCONDITION_HPP_

#include "FieldSpecificationBase.hpp"

namespace geos
{

/**
 * @class SourceFluxBoundaryCondition
 * A class to manage Neumann boundary conditions
 */
class SourceFluxBoundaryCondition : public FieldSpecificationBase
{
public:
  /**
   * @brief constructor
   * @param name the name of the FieldSpecificationBase in the data repository
   * @param parent the parent group of this group.
   */
  SourceFluxBoundaryCondition( string const & name, dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  SourceFluxBoundaryCondition() = delete;

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "SourceFlux"; }

  virtual const string getCatalogName() const override
  {
    return SourceFluxBoundaryCondition::catalogName();
  }

};

} /* namespace geos */

#endif /* GEOS_FIELDSPECIFICATION_SOURCEFLUXBOUNDARYCONDITION_HPP_ */
