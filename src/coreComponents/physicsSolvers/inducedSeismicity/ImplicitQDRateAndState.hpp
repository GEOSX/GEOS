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

#ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_IMPLICITQDRATEANDSTATE_HPP
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_IMPLICITQDRATEANDSTATE_HPP

#include "QDRateAndStateBase.hpp"

namespace geos
{

class ImplicitQDRateAndState : public QDRateAndStateBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  ImplicitQDRateAndState() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  ImplicitQDRateAndState( const string & name,
                          Group * const parent );

  /// Destructor
  virtual ~ImplicitQDRateAndState() override;

  static string derivedSolverPrefix() { return "Implicit";};

  struct viewKeyStruct : public QDRateAndStateBase::viewKeyStruct
  {
    /// target slip increment
    constexpr static char const * targetSlipIncrementString() { return "targetSlipIncrement"; }
  };

  virtual real64 setNextDt( real64 const & currentDt,
                            DomainPartition & domain ) override final;

  /**
   * @brief save the old state
   * @param subRegion
   */
  void updateSlip( ElementSubRegionBase & subRegion, real64 const dt ) const;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override final;
protected:

  void solveRateAndStateEquations( real64 const time_n,
                                   real64 const dt,
                                   DomainPartition & domain ) const;                               

  /// target slip rate
  real64 m_targetSlipIncrement;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_IMPLICITQDRATEANDSTATE_HPP */
