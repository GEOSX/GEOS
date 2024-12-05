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

/// THIS is an alternative implementation to avoid the use of the TractionUpdateWrapper

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEARTHQUAKE_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEARTHQUAKE_HPP

#include "physicsSolvers/inducedSeismicity/QuasiDynamicEQBase.hpp"

namespace geos
{

class QuasiDynamicEarthQuake : public QuasiDynamicEQBase 
{
public:

  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  QuasiDynamicEarthQuake() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  QuasiDynamicEarthQuake( const string & name,
                      Group * const parent );

  /// Destructor
  virtual ~QausiDynamicEarthQuake() override;

  static string catalogName() { return "QuasiDynamicEarthQuake"; }

  /**
   * @return Get the final class Catalog name
   */
  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public QuasiDynamicEQBase::viewKeyStruct
  {
    /// stress solver name
    static constexpr char const * contactSolverNameString() { return "stressSolverName"; }
  };

private:

PhysicsSolverBase * m_stressSolver;

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEQBASE_HPP */
