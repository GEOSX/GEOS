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

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_SPRINGSLIDER_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_SPRINGSLIDER_HPP

#include "physicsSolvers/QuasiDynamicEQBase.hpp"

namespace geos
{

class OneWayCoupledQuasiDynamicEQ : public QuasiDynamicEQBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  OneWayCoupledQuasiDynamicEQ() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  OneWayCoupledQuasiDynamicEQ( const string & name,
                               Group * const parent );

  /// Destructor
  virtual ~OneWayCoupledQuasiDynamicEQ() override;

  static string catalogName() { return "OneWayCoupledQuasiDynamicEQ"; }

  /**
   * @return Get the final class Catalog name
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override;

  struct viewKeyStruct : public PhysicsSolverBase::viewKeyStruct
  {
    /// Friction law name string
    constexpr static char const * frictionLawNameString() { return "frictionLawName"; }
  };

private:

  real64 updateStresses( real64 const & time_n,
                         real64 const & dt,
                         const int cycleNumber,
                         DomainPartition & domain ) const override final;

  virtual void postInputInitialization() override;


};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SPRINGSLIDER_HPP */
