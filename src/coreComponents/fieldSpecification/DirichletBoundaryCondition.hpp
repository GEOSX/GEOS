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
 * @file DirichletBoundaryCondition.hpp
 */

#ifndef GEOS_FIELDSPECIFICATION_DIRICHLETBOUNDARYCONDITION_HPP_
#define GEOS_FIELDSPECIFICATION_DIRICHLETBOUNDARYCONDITION_HPP_

#include "FieldSpecificationBase.hpp"

namespace geos
{

/**
 * @class DirichletBoundaryCondition
 * A class to manage Dirichlet boundary conditions
 */
class DirichletBoundaryCondition : public FieldSpecificationBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  DirichletBoundaryCondition( string const & name, dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  DirichletBoundaryCondition() = delete;

  /**
   * @brief destructor
   */
  virtual ~DirichletBoundaryCondition();

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "Dirichlet"; }

  virtual const string getCatalogName() const
  {
    return DirichletBoundaryCondition::catalogName();
  }

};

} /* namespace geos */

#endif /* GEOS_FIELDSPECIFICATION_DIRICHLETBOUNDARYCONDITION_HPP_ */
