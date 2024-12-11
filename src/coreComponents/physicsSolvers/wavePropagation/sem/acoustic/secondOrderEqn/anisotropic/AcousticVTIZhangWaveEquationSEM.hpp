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
 * @file AcousticVTIZhangWaveEquationSEM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIZHANGWAVEEQUATIONSEM_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIZHANGWAVEEQUATIONSEM_HPP_

#include "physicsSolvers/wavePropagation/shared/WaveSolverBase.hpp"
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticFields.hpp"
#include "AcousticVTIFields.hpp"

namespace geos
{

class AcousticVTIZhangWaveEquationSEM : public WaveSolverBase
{
public:

  using EXEC_POLICY = parallelDevicePolicy<  >;
  using ATOMIC_POLICY = AtomicPolicy< EXEC_POLICY >;

  AcousticVTIZhangWaveEquationSEM( const std::string & name,
                           Group * const parent );

  virtual ~AcousticVTIZhangWaveEquationSEM() override;

  AcousticVTIZhangWaveEquationSEM() = delete;
  AcousticVTIZhangWaveEquationSEM( AcousticVTIZhangWaveEquationSEM const & ) = delete;
  AcousticVTIZhangWaveEquationSEM( AcousticVTIZhangWaveEquationSEM && ) = default;

  AcousticVTIZhangWaveEquationSEM & operator=( AcousticVTIZhangWaveEquationSEM const & ) = delete;
  AcousticVTIZhangWaveEquationSEM & operator=( AcousticVTIZhangWaveEquationSEM && ) = delete;

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "acoustic"; }

  static string catalogName() { return "AcousticVTIZhangSEM"; }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

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

  virtual real64 explicitStepBackward( real64 const & time_n,
                                       real64 const & dt,
                                       integer const cycleNumber,
                                       DomainPartition & domain,
                                       bool const computeGradient ) override;

  /**@}*/

  /**
   * @brief Multiply the precomputed term by the Ricker and add to the right-hand side
   * @param cycleNumber the cycle number/step number of evaluation of the source
   * @param rhs the right hand side vector to be computed
   */
  virtual void addSourceToRightHandSide( real64 const & time_n, arrayView1d< real32 > const rhs );

  void precomputeSurfaceFieldIndicator( DomainPartition & domain );


  /**
   * @brief Initialize Perfectly Matched Layer (PML) information
   */
  virtual void initializePML() override;

  /**
   */
  virtual real64 computeTimeStep( real64 & dtOut ) override;



  /**
   * @brief Overridden from ExecutableGroup. Used to write last seismogram if needed.
   */
  virtual void cleanup( real64 const time_n, integer const cycleNumber, integer const eventCounter, real64 const eventProgress, DomainPartition & domain ) override;

  struct viewKeyStruct : WaveSolverBase::viewKeyStruct
  {
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
                               DomainPartition & domain,
                               bool const isForward );

  void computeUnknowns( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain,
                        MeshLevel & mesh,
                        arrayView1d< string const > const & regionNames,
                        bool const isForward );

  void synchronizeUnknowns( real64 const & time_n,
                            real64 const & dt,
                            DomainPartition & domain,
                            MeshLevel & mesh,
                            arrayView1d< string const > const & regionNames );

  void prepareNextTimestep( MeshLevel & mesh );

protected:

  virtual void postInputInitialization() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  /**
   * @brief Locate sources and receivers position in the mesh elements, evaluate the basis functions at each point and save them to the
   * corresponding elements nodes.
   * @param baseMesh the level-0 mesh
   * @param mesh mesh of the computational domain
   */
  virtual void precomputeSourceAndReceiverTerm( MeshLevel & baseMesh, MeshLevel & mesh, arrayView1d< string const > const & regionNames ) override;

  /**
   * @brief Apply free surface condition to the face define in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyFreeSurfaceBC( real64 const time, DomainPartition & domain ) override;

  /**
   * @brief Apply Perfectly Matched Layer (PML) to the regions defined in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyPML( real64 const time, DomainPartition & domain ) override;

  /// Pressure_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_pressureNp1AtReceivers;

  /// Array of size the number of receivers and filled with 0.5 (used for calculating the seismos)
  array1d< real32 > m_seismoCoeff;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIZHANGWAVEEQUATIONSEM_HPP_ */
