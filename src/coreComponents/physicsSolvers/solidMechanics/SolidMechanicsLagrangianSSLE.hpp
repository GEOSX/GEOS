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
 * @file SolidMechanicsLagrangianSSLE.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_

#include "SolidMechanicsLagrangianFEM.hpp"

namespace geos
{

/**
 * @class SolidMechanicsLagrangianSSLE
 *
 * This class contains an implementation of a small strain linear elastic solution to the equations of motion which are
 * called through the interface in SolidMechanicsLagrangianFEM.
 */
class SolidMechanicsLagrangianSSLE : public SolidMechanicsLagrangianFEM
{
public:
  SolidMechanicsLagrangianSSLE( string const & name,
                                Group * const parent );
  virtual ~SolidMechanicsLagrangianSSLE() override;

  static string catalogName() { return "SolidMechanicsLagrangianSSLE"; }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }


};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_ */
