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

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEQRK32_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEQRK32_HPP

#include "physicsSolvers/SolverBase.hpp"
#include "kernels/RateAndStateKernels.hpp"

namespace geos
{

class QuasiDynamicEQRK32 : public SolverBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  QuasiDynamicEQRK32() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  QuasiDynamicEQRK32( const string & name,
                  Group * const parent );

  /// Destructor
  virtual ~QuasiDynamicEQRK32() override;

  static string catalogName() { return "QuasiDynamicEQRK32"; }

  /**
   * @return Get the final class Catalog name
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override;

  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    /// stress solver name
    static constexpr char const * stressSolverNameString() { return "stressSolverName"; }
    /// Friction law name string
    constexpr static char const * frictionLawNameString() { return "frictionLawName"; }
    /// Friction law name string
    constexpr static char const * shearImpedanceString() { return "shearImpedance"; }
    /// target slip increment
    constexpr static char const * timeStepTol() { return "timeStepTol"; }
  };

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override final;

private:

  virtual real64 setNextDt( real64 const & currentDt,
                            DomainPartition & domain ) override final;

  virtual void postInputInitialization() override;

  void stepRateStateODESubstage( integer const stageIndex,
                                                     real64 const dt,
                                                     DomainPartition & domain ) const;

  void stepRateStateODEAndComputeError(real64 const dt,
                                                             DomainPartition & domain ) const;

  real64 updateStresses( real64 const & time_n,
                         real64 const & dt,
                         const int cycleNumber,
                         DomainPartition & domain ) const;

   /**
   * @brief updates rate-and-state slip velocity and slip rate
   * @param domain
   */
  void updateSlipVelocity( real64 const & time_n,
                           real64 const & dt,
                           DomainPartition & domain ) const;                        

  /**
   * @brief save the current state
   * @param domain
   */
  void saveState( DomainPartition & domain) const;


  /// pointer to stress solver
  SolverBase * m_stressSolver;

  /// stress solver name
  string m_stressSolverName;

  /// shear impedance
  real64 m_shearImpedance;

  /// time step tolerance
  real64 m_timestepAbsTol;

  real64 m_timestepRelTol;

  real64 m_timestepAcceptSafety;
  
  real64 m_prevTimestepErrors[2];

  /// time step error
  real64 m_beta[3];

  bool m_successfulStep;

  rateAndStateKernels::RK32Table m_butcherTable;

  class SpringSliderParameters
  {
public:

    GEOS_HOST_DEVICE
    SpringSliderParameters( real64 const normalTraction, real64 const a, real64 const b, real64 const Dc ):
      tauRate( 1e-4 ),
      springStiffness( 0.0 )
    {
      real64 const criticalStiffness = normalTraction * (b - a) / Dc;
      springStiffness = 0.9 * criticalStiffness;
    }

    /// Default copy constructor
    SpringSliderParameters( SpringSliderParameters const & ) = default;

    /// Default move constructor
    SpringSliderParameters( SpringSliderParameters && ) = default;

    /// Deleted default constructor
    SpringSliderParameters() = delete;

    /// Deleted copy assignment operator
    SpringSliderParameters & operator=( SpringSliderParameters const & ) = delete;

    /// Deleted move assignment operator
    SpringSliderParameters & operator=( SpringSliderParameters && ) =  delete;

    real64 tauRate;

    real64 springStiffness;
  };
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEQRK32_HPP */
