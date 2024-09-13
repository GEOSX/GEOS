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

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_SEISMICITYRATE_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_SEISMICITYRATE_HPP

#include "codingUtilities/EnumStrings.hpp"   // facilities for enum-string conversion (for reading enum values from XML input)
#include "physicsSolvers/SolverBase.hpp"  // an abstraction class shared by all physics solvers
#include "fieldSpecification/FieldSpecificationManager.hpp" // a manager that can access and set values on the discretized domain

#include "physicsSolvers/inducedSeismicity/inducedSeismicityFields.hpp"

#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

namespace geos
{

class SeismicityRate : public SolverBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  SeismicityRate() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  SeismicityRate( const string & name,
                  Group * const parent );

  /// Destructor
  virtual ~SeismicityRate() override;

  static string catalogName() { return "SeismicityRate"; }

  /**
   * @return Get the final class Catalog name
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override;

  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    static constexpr char const * stressSolverNameString() { return "stressSolverName"; }
    static constexpr char const * initialFaultNormalTractionString() { return "initialFaultNormalTraction"; }
    static constexpr char const * initialFaultShearTractionString() { return "initialFaultShearTraction"; }
    static constexpr char const * faultNormalDirectionString() { return "faultNormalDirection"; }
    static constexpr char const * faultShearDirectionString() { return "faultShearDirection"; }
    static constexpr char const * directEffectString() { return "directEffect"; }
    static constexpr char const * backgroundStressingRateString() { return "backgroundStressingRate"; }
  };

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override final;

  /**
   * @brief single step advance in computing the seismicity rate based on
   *  stress history according to closed form integral solution (Heimisson & Segall, 2018)
   *  to the ODE formulated by Dieterich, 1994
   * @param time_n time at previous converged step
   * @param dt time step size
   * @param subRegion ElementSubRegionBase to compute the solution in
   */
  void integralSolverStep( real64 const & time_n,
                           real64 const & dt,
                           ElementSubRegionBase & subRegion );

  /**
   * @brief called in SolverStep after member stress solver is called to
   *  project the stress state to pre-defined fault orientations
   * @param subRegion The ElementSubRegionBase that will have the stress information
   */
  void updateFaultTraction( ElementSubRegionBase & subRegion ) const;

  /**
   * @brief save the old state
   * @param subRegion
   */
  void saveOldState( ElementSubRegionBase & subRegion ) const;

  /**
   * @brief
   * @param biotCoefficient
   * @param pres
   * @param sig
   * @param tau
   */
  void computeTotalStressOnFault( arrayView1d< real64 const > const biotCoefficient,
                                  arrayView1d< real64 const > const pres,
                                  real64 const (&faultNormalProjectionTensor)[6],
                                  real64 const (&faultShearProjectionTensor)[6],
                                  arrayView1d< real64 > const sig,
                                  arrayView1d< real64 > const tau ) const;

  /**
   * @brief called in SolverStep before member stress solver is called to
   *  project the initial stress state to pre-defined fault orientations,
   *  only when cycleNumber == 0
   * @param cycleNumber current cycle number
   * @param domain The DomainPartition of the problem
   */
  void initializeFaultTraction( real64 const time_n,
                                integer const cycleNumber,
                                DomainPartition & domain ) const;


protected:

  /**
   * @brief update the stresses either by asking the stressSolver or by applying b.c.
   * @param time_n the current time
   * @param dt the current time step
   * @param cycleNumber the cycle number
   * @param domain the DomainPartion group
   */
  real64 updateStresses( real64 const & time_n,
                         real64 const & dt,
                         const int cycleNumber,
                         DomainPartition & domain ) const;



  void constructFaultStressProjectionTensors( real64 ( &faultNormalProjectionTensor )[6],
                                              real64 ( &faultShearProjectionTensor )[6] ) const;

  virtual void postInputInitialization() override;

  /// pointer to stress solver
  SolverBase * m_stressSolver;

  /// stress solver name string
  string m_stressSolverName;

  /// fault orientation: normal direction
  R1Tensor m_faultNormalDirection;

  /// fault orientation: shear direction
  R1Tensor m_faultShearDirection;

  /// direct effect coefficient
  real64 m_directEffect;

  /// bacground stressing rate
  real64 m_backgroundStressingRate;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_SEISMICITYRATE_HPP */
