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
 * @file AcousticPOD.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICPOD_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICPOD_HPP_

#include "physicsSolvers/wavePropagation/shared/WaveSolverBase.hpp"
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticFields.hpp"

namespace geos
{

class AcousticPOD : public WaveSolverBase
{
public:

  using EXEC_POLICY = parallelDevicePolicy<  >;
  using ATOMIC_POLICY = AtomicPolicy< EXEC_POLICY >;

  AcousticPOD( const std::string & name,
               Group * const parent );

  virtual ~AcousticPOD() override;

  AcousticPOD() = delete;
  AcousticPOD( AcousticPOD const & ) = delete;
  AcousticPOD( AcousticPOD && ) = default;

  AcousticPOD & operator=( AcousticPOD const & ) = delete;
  AcousticPOD & operator=( AcousticPOD && ) = delete;

  static string catalogName() { return "AcousticPOD"; }

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
                                       bool const computeGradient ) override {return 0;};


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
                                       arrayView2d< real32 > varAtReceivers,
                                       arrayView1d< real32 > coeffs = {},
                                       bool add = false );

  void computeMassAndDampingPOD( arrayView2d< real32 > const massPOD,
				 arrayView2d< real32 > const massGradientPOD,
				 arrayView2d< real32 > const dampingPOD,
				 arrayView1d< real32 const > const mass,
				 arrayView1d< real32 const > const massGradient,
				 arrayView1d< real32 const > const damping,
				 int const countPhi,
				 int const shotIndex,
				 arrayView1d< localIndex const > const nodesGhostRank );


  void computeSeismoTracePOD( real64 const time_n,
			      real64 const dt,
			      real64 const timeSeismo,
			      localIndex iSeismo,
			      arrayView2d< real64 const > const receiverConstants,
			      arrayView1d< localIndex const > const receiverIsLocal,
			      localIndex const nsamplesSeismoTrace,
			      localIndex const outputSeismoTrace,
			      arrayView1d< real32 const > const var_np1,
			      arrayView1d< real32 const > const var_n,
			      arrayView2d< real32 > varAtReceivers );


  void computeInitialConditions();
  /**
   * @brief Initialize Perfectly Matched Layer (PML) information
   */
  virtual void initializePML() override;


  virtual real64 computeTimeStep( real64 & dtOut ) override;
  /**
   * @brief Overridden from ExecutableGroup. Used to write last seismogram if needed.
   */
  virtual void cleanup( real64 const time_n, integer const cycleNumber, integer const eventCounter, real64 const eventProgress, DomainPartition & domain ) override;

  struct viewKeyStruct : WaveSolverBase::viewKeyStruct
  {

    static constexpr char const * pressureNp1AtReceiversString() { return "pressureNp1AtReceivers"; }
    static constexpr char const * computePODmatrixString() { return "computePODmatrix"; }
    static constexpr char const * computeSourceValueString() { return "computeSourceValue"; }
    static constexpr char const * massPODString() { return "massPOD"; }
    static constexpr char const * massGradientPODString() { return "massGradientPOD"; }
    static constexpr char const * dampingPODString() { return "dampingPOD"; }
    static constexpr char const * invAPODString() { return "invAPOD"; }
    static constexpr char const * sourceConstantsPODString() { return "sourceConstantsPOD"; }
    static constexpr char const * a_np1String() { return "a_np1"; }
    static constexpr char const * a_nString() { return "a_n"; }
    static constexpr char const * a_nm1String() { return "a_nm1"; }
    static constexpr char const * perturbationString() { return "perturbation"; }
    static constexpr char const * sizePODString() { return "sizePOD"; }
    static constexpr char const * orderInitString() { return "orderInit"; }

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

protected:

  virtual void postInputInitialization() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  /**
   * @brief Locate sources and receivers position in the mesh elements, evaluate the basis functions at each point and save them to the
   * corresponding elements nodes.
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

  localIndex getNumNodesPerElem();

  /// Pressure_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_pressureNp1AtReceivers;

  /// Taylor order to compute the initial conditions
  localIndex m_orderInit;

  /// Whether or not to compute the POD matrices
  int m_computePODmatrix;

  /// Whether or not to compute the source and receivers
  int m_computeSourceValue;

  ///Mass Matrix in the POD basis forward
  array2d< real32 > m_massPOD;

  ///Mass Matrix with gradient in the pod basis
  array2d< real32 > m_massGradientPOD;

  ///Damping Matrix in the POD basis forward
  array2d< real32 > m_dampingPOD;

  ///Inverse scalar product POD Matrix
  array2d< real32 > m_invAPOD;

  // Model perturbation coefficient
  real32 m_perturbation;

  // Coefficient of the solution in the POD basis at time n+1
  array1d< real32 > m_a_np1;

  // Coefficient of the solution in the POD basis at time n
  array1d< real32 > m_a_n;

  // Coefficient of the solution in the POD basis at time n-1
  array1d< real32 > m_a_nm1;

  // Second derivative of a for all dtWaveField time steps
  array2d< real32 > m_a_dt2;

  // Right hand side in the POD basis
  array1d< real32 > m_rhsPOD;

  // Contribution on POD basis for the source
  array2d< real64 > m_sourceConstantsPOD;

  // Contribution on POD basis for the receivers
  array2d< real64 > m_receiverConstantsPOD;

  int m_sizePOD;
};


namespace fields
{

DECLARE_FIELD( Pressure_nm1,
               "pressure_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n-1." );

DECLARE_FIELD( Pressure_n,
               "pressure_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n." );

DECLARE_FIELD( Pressure_np1,
               "pressure_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar pressure at time n+1." );

DECLARE_FIELD( ForcingRHS,
               "rhs",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS" );

DECLARE_FIELD( AcousticMassVector,
               "acousticMassVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Mass Matrix." );

DECLARE_FIELD( AcousticMassGradientVector,
               "acousticMassGradientVector",
	       array1d< real32 >,
               0,
               NOPLOT,
	       WRITE_AND_READ,
	       "Diagonal of the Mass Matrix with gradient of velocity." );

DECLARE_FIELD( AcousticDampingVector,
               "acousticDampingVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix." );

DECLARE_FIELD( AcousticVelocity,
               "acousticVelocity",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Medium velocity of the cell" );

DECLARE_FIELD( AcousticDensity,
               "acousticDensity",
	       array1d< real32 >,
	       0,
               NOPLOT,
               WRITE_AND_READ,
	       "Medium density of the cell" );

DECLARE_FIELD( AcousticFreeSurfaceFaceIndicator,
               "acousticFreeSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
	       WRITE_AND_READ,
               "Free surface indicator, 1 if a face is on free surface 0 otherwise." );

DECLARE_FIELD( AcousticFreeSurfaceNodeIndicator,
               "acousticFreeSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
	       WRITE_AND_READ,
               "Free surface indicator, 1 if a node is on free surface 0 otherwise." );

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

DECLARE_FIELD( AuxiliaryVar1PML,
               "auxiliaryVar1PML",
               array2d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "PML vectorial auxiliary variable 1." );

DECLARE_FIELD( AuxiliaryVar2PML,
               "auxiliaryVar2PML",
               array2d< real32 >,
               0,
               NOPLOT,
               NO_WRITE,
               "PML vectorial auxiliary variable 2." );

DECLARE_FIELD( AuxiliaryVar3PML,
               "auxiliaryVar3PML",
               array1d< real32 >,
               0,
               NOPLOT,
               NO_WRITE,
               "PML scalar auxiliary variable 3." );

DECLARE_FIELD( AuxiliaryVar4PML,
               "auxiliaryVar4PML",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "PML scalar auxiliary variable 4." );

}


} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICPOD_HPP_ */
