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
 * @file AcousticFirstOrderWaveEquationSEM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIRSTORDERWAVEEQUATIONSEM_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIRSTORDERWAVEEQUATIONSEM_HPP_

#include "mesh/MeshFields.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticFields.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverBase.hpp"

namespace geos
{

class AcousticFirstOrderWaveEquationSEM : public WaveSolverBase
{
public:

  using EXEC_POLICY = parallelDevicePolicy< >;
  using ATOMIC_POLICY = parallelDeviceAtomic;

  AcousticFirstOrderWaveEquationSEM( const std::string & name,
                                     Group * const parent );

  virtual ~AcousticFirstOrderWaveEquationSEM() override;

  static string catalogName() { return "AcousticFirstOrderSEM"; }
  /**
   * @copydoc SolverBase::getCatalogName()
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

  /**
   * @brief Method to get the maximum wavespeed on a mesh (usually the P-wavespeed)
   */
  virtual real32 getGlobalMaxWavespeed( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) override;

  /**
   * @brief Initialize Perfectly Matched Layer (PML) information
   */
  virtual void initializePML() override;


  /**
   * @brief Overridden from ExecutableGroup. Used to write last seismogram if needed.
   */
  virtual void cleanup( real64 const time_n, integer const cycleNumber, integer const eventCounter, real64 const eventProgress, DomainPartition & domain ) override;


  struct viewKeyStruct : WaveSolverBase::viewKeyStruct
  {

    static constexpr char const * pressureNp1AtReceiversString() { return "pressureNp1AtReceivers"; }

    static constexpr char const * uxNp1AtReceiversString() { return "uxNp1AtReceivers"; }
    static constexpr char const * uyNp1AtReceiversString() { return "uyNp1AtReceivers"; }
    static constexpr char const * uzNp1AtReceiversString() { return "uzNp1AtReceivers"; }

    static constexpr char const * sourceElemString() { return "sourceElem"; }
    static constexpr char const * sourceRegionString() { return "sourceRegion"; }

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

  virtual real64 computeTimeStep( real64 & dtOut ) override;

  /// Pressure_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_pressureNp1AtReceivers;

  /// Pressure_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_uxNp1AtReceivers;

  /// Pressure_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_uyNp1AtReceivers;

  /// Pressure_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_uzNp1AtReceivers;

  /// Array containing the elements which contain a source
  array1d< localIndex > m_sourceElem;

  /// Array containing the elements which contain the region which the source belongs
  array1d< localIndex > m_sourceRegion;
};


} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIRSTORDERWAVEEQUATIONSEM_HPP_ */
