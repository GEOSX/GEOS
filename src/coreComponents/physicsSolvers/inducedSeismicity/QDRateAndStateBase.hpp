/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
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
   * @brief Save the current state of the solver fields
   * @param domain the domain object
   */
  void saveState( DomainPartition & domain ) const;

  /**
   * @brief Check that only one of slip rate or slip velocity are specified as initial conditions
   * and initialize the unspecified field
   * @param subRegion the element subregion
   */
  void enforceRateAndVelocityConsistency( SurfaceElementSubRegion & subRegion ) const;

  /**
   * @brief Compute stresses and update tractions on the fault
   * @param time_n the current time
   * @param dt the time step
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   */
  virtual real64 updateStresses( real64 const & time_n,
                                 real64 const & dt,
                                 const int cycleNumber,
                                 DomainPartition & domain ) const = 0;

  /**
   * @brief Apply initial conditions to fields on the fault
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   */
  virtual void applyInitialConditionsToFault( int const cycleNumber,
                                              DomainPartition & domain ) const;

protected:

  /// shear impedance
  real64 m_shearImpedance;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCEDSEISMICTY_QDRATEANDSTATEBASE_HPP */
