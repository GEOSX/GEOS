/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file ElasticWaveEquationSEM.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEM_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEM_HPP_

#include "physicsSolvers/wavePropagation/shared/WaveSolverBase.hpp"
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "physicsSolvers/wavePropagation/sem/elastic/shared/ElasticFields.hpp"

namespace geos
{

class ElasticWaveEquationSEM : public WaveSolverBase
{
public:

  using EXEC_POLICY = parallelDevicePolicy<  >;
  using ATOMIC_POLICY = parallelDeviceAtomic;

  ElasticWaveEquationSEM( const std::string & name,
                          Group * const parent );

  virtual ~ElasticWaveEquationSEM() override;

  ElasticWaveEquationSEM() = delete;
  ElasticWaveEquationSEM( ElasticWaveEquationSEM const & ) = delete;
  ElasticWaveEquationSEM( ElasticWaveEquationSEM && ) = default;

  ElasticWaveEquationSEM & operator=( ElasticWaveEquationSEM const & ) = delete;
  ElasticWaveEquationSEM & operator=( ElasticWaveEquationSEM && ) = delete;

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "elastic"; }

  static string catalogName() { return "ElasticSEM"; }
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
  /**@}*/

  /**
   * @brief Multiply the precomputed term by the Ricker and add to the right-hand side
   * @param cycleNumber the cycle number/step number of evaluation of the source
   * @param rhsx the right hand side vector to be computed (x-component)
   * @param rhsy the right hand side vector to be computed (y-component)
   * @param rhsz the right hand side vector to be computed (z-component)
   */
  void addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhsx, arrayView1d< real32 > const rhsy, arrayView1d< real32 > const rhsz );

  /**
   * TODO: move implementation into WaveSolverBase once 'm_receiverIsLocal' is also moved
   * @brief Compute DAS data as a difference of the field at two points, from the appropriate three-component receiver pairs, when the DAS
   * type is set to 2
   * @param xCompRcv the array holding the x-component of pairs of receivers
   * @param yCompRcv the array holding the y-component of pairs of receivers
   * @param zCompRcv the array holding the z-component of pairs of receivers
   */
  void computeDAS( arrayView2d< real32 > const xCompRcv,
                   arrayView2d< real32 > const yCompRcv,
                   arrayView2d< real32 > const zCompRcv );

  /**
   * @brief Overridden from ExecutableGroup. Used to write last seismogram if needed.
   */
  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  struct viewKeyStruct : WaveSolverBase::viewKeyStruct
  {
    static constexpr char const * displacementXNp1AtReceiversString() { return "displacementXNp1AtReceivers"; }
    static constexpr char const * displacementYNp1AtReceiversString() { return "displacementYNp1AtReceivers"; }
    static constexpr char const * displacementZNp1AtReceiversString() { return "displacementZNp1AtReceivers"; }

    static constexpr char const * dasSignalNp1AtReceiversString() { return "dasSignalNp1AtReceivers"; }

    static constexpr char const * sourceForceString() { return "sourceForce"; }
    static constexpr char const * sourceMomentString() { return "sourceMoment"; }

    static constexpr char const * useVtiString() { return "useVTI"; }

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

  void computeUnknowns( real64 const & time_n,
                        real64 const & dt,
                        integer const cycleNumber,
                        DomainPartition & domain,
                        MeshLevel & mesh,
                        arrayView1d< string const > const & regionNames );

  void synchronizeUnknowns( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain,
                            MeshLevel & mesh,
                            arrayView1d< string const > const & regionNames );

  void prepareNextTimestep( MeshLevel & mesh );

  /**
   * @brief Computes the minimum attenuation quality factor over all the mesh. This is useful for computing anelasticity coefficients, which
   * are usually global parameters
   */
  real32 computeGlobalMinQFactor();

protected:

  virtual void postInputInitialization() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  /**
   * @brief Locate sources and receivers position in the mesh elements, evaluate the basis functions at each point and save them to the
   * corresponding elements nodes.
   * @param mesh mesh of the computational domain
   * @param regionNames the names of the region you loop on
   */
  virtual void precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) override;

  /**
   * @brief Apply free surface condition to the face define in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyFreeSurfaceBC( real64 const time, DomainPartition & domain ) override;

  /**
   * @brief Initialize Perfectly Matched Layer (PML) information
   */
  virtual void initializePML() override;

  /**
   * @brief Apply Perfectly Matched Layer (PML) to the regions defined in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyPML( real64 const time, DomainPartition & domain ) override;

  /// Constant part of the source for the nodes listed in m_sourceNodeIds in x-direction
  array2d< real64 > m_sourceConstantsx;

  /// Constant part of the source for the nodes listed in m_sourceNodeIds in y-direction
  array2d< real64 > m_sourceConstantsy;

  /// Constant part of the source for the nodes listed in m_sourceNodeIds in z-direction
  array2d< real64 > m_sourceConstantsz;

  /// Displacement_np1 at the receiver location for each time step for each receiver (x-component)
  array2d< real32 > m_displacementXNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver (y-component)
  array2d< real32 > m_displacementYNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver (z-component)
  array2d< real32 > m_displacementZNp1AtReceivers;

  /// DAS receiver signal at np1 for each time step for each receiver (z-component)
  array2d< real32 > m_dasSignalNp1AtReceivers;

  /// Vector describing the force of the source
  R1Tensor m_sourceForce;

  /// Symmetric tensor describing the moment of the source
  R2SymTensor m_sourceMoment;

  /// Flag to appliy VTI anisotropy
  integer m_useVTI;

};

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASSTICWAVEEQUATIONSEM_HPP_ */
