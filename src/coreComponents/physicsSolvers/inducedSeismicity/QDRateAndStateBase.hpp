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

#ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_QDRATEANDSTATEBASE_HPP
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_QDRATEANDSTATEBASE_HPP

#include "physicsSolvers/PhysicsSolverBase.hpp"

namespace geos
{

class QDRateAndStateBase : public PhysicsSolverBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  QDRateAndStateBase() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  QDRateAndStateBase( const string & name,
                      Group * const parent );

  /// Destructor
  virtual ~QDRateAndStateBase() override;

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override;

  struct viewKeyStruct : public PhysicsSolverBase::viewKeyStruct
  {
    /// Friction law name string
    constexpr static char const * frictionLawNameString() { return "frictionLawName"; }
    /// Friction law name string
    constexpr static char const * shearImpedanceString() { return "shearImpedance"; }
  };

  /**
   * @brief save the current state
   * @param domain
   */
  void saveState( DomainPartition & domain ) const;

protected:

  virtual real64 updateStresses( real64 const & time_n,
                                 real64 const & dt,
                                 const int cycleNumber,
                                 DomainPartition & domain ) const = 0;

  virtual void applyInitialConditionsToFault( int const cycleNumber,
                                              DomainPartition & domain ) const;

  /// shear impedance
  real64 m_shearImpedance;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCEDSEISMICTY_QDRATEANDSTATEBASE_HPP */
