/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOS Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file AcousticVTIWaveEquationSEM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_AcousticVTIWaveEquationSEM_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_AcousticVTIWaveEquationSEM_HPP_

#include "WaveSolverBase.hpp"
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/SolverBase.hpp"

namespace geos
{

class AcousticVTIWaveEquationSEM : public WaveSolverBase
{
public:

  using EXEC_POLICY = parallelDevicePolicy< 32 >;
  using ATOMIC_POLICY = AtomicPolicy< EXEC_POLICY >;

  AcousticVTIWaveEquationSEM( const std::string & name,
                           Group * const parent );

  virtual ~AcousticVTIWaveEquationSEM() override;

  AcousticVTIWaveEquationSEM() = delete;
  AcousticVTIWaveEquationSEM( AcousticVTIWaveEquationSEM const & ) = delete;
  AcousticVTIWaveEquationSEM( AcousticVTIWaveEquationSEM && ) = default;

  AcousticVTIWaveEquationSEM & operator=( AcousticVTIWaveEquationSEM const & ) = delete;
  AcousticVTIWaveEquationSEM & operator=( AcousticVTIWaveEquationSEM && ) = delete;


  static string catalogName() { return "AcousticSEM_VTI"; }

  bool dirExists( const std::string & dirName );

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & meshBodies ) override final;


  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual real64 explicitStepForward( real64 const & time_n,
                                      real64 const & dt,
                                      integer const cycleNumber,
                                      DomainPartition & domain,
                                      bool const computeGradient ) override;


  virtual real64 explicitStepBackward( real64 const & GEOS_UNUSED_PARAM(time_n),
                                       real64 const & GEOS_UNUSED_PARAM(dt),
                                       integer const GEOS_UNUSED_PARAM(cycleNumber),
                                       DomainPartition & GEOS_UNUSED_PARAM(domain),
                                       bool const GEOS_UNUSED_PARAM(computeGradient) ) override;

  /**@}*/

  /**
   * @brief Multiply the precomputed term by the Ricker and add to the right-hand side
   * @param cycleNumber the cycle number/step number of evaluation of the source
   * @param rhs the right hand side vector to be computed
   */
  virtual void addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhs );

  /**
   * TODO: move implementation into WaveSolverBase
   * @brief Computes the traces on all receivers (see @computeSeismoTraces) up to time_n+dt
   * @param time_n the time corresponding to the field values pressure_n
   * @param dt the simulation timestep
   * @param var_np1 the field values at time_n + dt
   * @param var_n the field values at time_n
   * @param varAtReceivers the array holding the trace values, where the output is written
   */
  virtual void computeAllSeismoTraces( real64 const time_n,
                                       real64 const dt,
                                       arrayView1d< real32 const > const var_np1,
                                       arrayView1d< real32 const > const var_n,
                                       arrayView2d< real32 > varAtReceivers );

  /**
   * @brief Overridden from ExecutableGroup. Used to write last seismogram if needed.
   */
  virtual void cleanup( real64 const time_n, integer const cycleNumber, integer const eventCounter, real64 const eventProgress, DomainPartition & domain ) override;

  struct viewKeyStruct : WaveSolverBase::viewKeyStruct
  {
    static constexpr char const * sourceNodeIdsString() { return "sourceNodeIds"; }
    static constexpr char const * sourceConstantsString() { return "sourceConstants"; }
    static constexpr char const * sourceIsAccessibleString() { return "sourceIsAccessible"; }

    static constexpr char const * receiverNodeIdsString() { return "receiverNodeIds"; }
    static constexpr char const * receiverConstantsString() {return "receiverConstants"; }
    static constexpr char const * receiverIsLocalString() { return "receiverIsLocal"; }

    static constexpr char const * pressureNp1AtReceiversString() { return "pressureNp1AtReceivers"; }

  } waveEquationViewKeys;


  /** internal function to the class to compute explicitStep either for backward or forward.
   * (requires not to be private because it is called from GEOS_HOST_DEVICE method)
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @return return the timestep that was achieved during the step.
   */
  real64 explicitStepInternal( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain );

  /**
   * @brief (Empty but must be defined) Initialize Perfectly Matched Layer (PML) information
   */
  virtual void initializePML() override {}

  /**
   * @brief  (Empty but must be defined) Apply Perfectly Matched Layer (PML) to the regions defined in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyPML( real64 const GEOS_UNUSED_PARAM(time), DomainPartition & GEOS_UNUSED_PARAM(domain) ) override {}

protected:

  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  /**
   * @brief Locate sources and receivers position in the mesh elements, evaluate the basis functions at each point and save them to the
   * corresponding elements nodes.
   * @param mesh mesh of the computational domain
   */
  virtual void precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) override;

  /**
   * @brief Compute the surface Field indicators from the xml
   * @param time the time to apply the BC TODO: a priori useless...
   * @param domain the partition domain
   */
  virtual void precomputeSurfaceFieldIndicator( real64 const time, DomainPartition & domain );

  /**
   * @brief Apply free surface condition to the face define in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyFreeSurfaceBC( real64 const time, DomainPartition & domain ) override;

  /// Pressure_p_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_pressureNp1AtReceivers;

};


namespace fields
{

DECLARE_FIELD( Pressure_p_nm1,
               "pressure_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n-1." );

DECLARE_FIELD( Pressure_p_n,
               "pressure_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n." );

DECLARE_FIELD( Pressure_p_np1,
               "pressure_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar pressure at time n+1." );

DECLARE_FIELD( Pressure_q_nm1,
               "pressure_q_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar auxiliary pressure q at time n-1." );

DECLARE_FIELD( Pressure_q_n,
               "pressure_q_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar auxiliary pressure q at time n." );

DECLARE_FIELD( Pressure_q_np1,
               "pressure_q_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar auxiliary pressure q at time n+1." );

DECLARE_FIELD( ForcingRHS,
               "rhs",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS" );

DECLARE_FIELD( MassVector,
               "massVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Mass Matrix." );

DECLARE_FIELD( DampingVector_p,
               "dampingVector_p",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix for p terms in p equation." );

DECLARE_FIELD( DampingVector_pq,
               "dampingVector_pq",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix for q terms in p equation." );

DECLARE_FIELD( DampingVector_q,
               "dampingVector_q",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix for q terms in q equation." );

DECLARE_FIELD( DampingVector_qp,
               "dampingVector_qp",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix for p terms in q equation." );

DECLARE_FIELD( MediumVelocity,
               "mediumVelocity",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Medium velocity of the cell" );

DECLARE_FIELD( Delta,
               "delta",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Delta thomsen anisotropy parameter" );

DECLARE_FIELD( Epsilon,
               "epsilon",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Epsilon thomsen anisotropy parameter" );

DECLARE_FIELD( F,
               "f",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "f quantity in VTI/TTI Fletcher's equations" );

DECLARE_FIELD( StiffnessVector_p,
               "stiffnessVector_p",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Stiffness vector contains R_h*Pressure_n." );

DECLARE_FIELD( StiffnessVector_q,
               "stiffnessVector_q",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Stiffness vector contains R_h*Pressure_n." );

DECLARE_FIELD( FreeSurfaceFaceIndicator,
               "freeSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a face is on free surface 0 otherwise." );

DECLARE_FIELD( FreeSurfaceNodeIndicator,
               "freeSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a node is on free surface 0 otherwise." );

DECLARE_FIELD( LateralSurfaceFaceIndicator,
               "lateralSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a face is on a lateral surface 0 otherwise." );

DECLARE_FIELD( LateralSurfaceNodeIndicator,
               "lateralSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Lateral surface indicator, 1 if a face is on a lateral surface 0 otherwise." );

DECLARE_FIELD( BottomSurfaceFaceIndicator,
               "bottomSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Bottom surface indicator, 1 if a face is on the bottom surface 0 otherwise." );

DECLARE_FIELD( BottomSurfaceNodeIndicator,
               "bottomSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Bottom surface indicator, 1 if a face is on the bottom surface 0 otherwise." );

DECLARE_FIELD( PressureDoubleDerivative,
               "pressureDoubleDerivative",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Double derivative of the pressure for each node to compute the gradient" );

DECLARE_FIELD( PartialGradient,
               "partialGradient",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Partiel gradient computed during backward propagation" );

//Debug


DECLARE_FIELD( StiffnessVector,
               "stiffnessVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Stiffness vector contains R_h*Pressure_n." );
               
DECLARE_FIELD( DampingVector,
               "dampingVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix." );
}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_AcousticVTIWaveEquationSEM_HPP_ */
