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
 * @file SolidMechanicsMPM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "common/TimingMacros.hpp"
#include "kernels/SolidMechanicsLagrangianFEMKernels.hpp"
#include "kernels/ExplicitMPM.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "MPMSolverFields.hpp"

namespace geos
{

class SpatialPartition;


/**
 * @class SolidMechanicsMPM
 *
 * This class implements a material point method solution to the equations of motion.
 */
class SolidMechanicsMPM : public SolverBase
{
public:

  /**
   * @enum TimeIntegrationOption
   *
   * The options for time integration
   */
  enum class TimeIntegrationOption : integer
  {
    QuasiStatic,      //!< QuasiStatic
    ImplicitDynamic,  //!< ImplicitDynamic
    ExplicitDynamic   //!< ExplicitDynamic
  };

  /**
   * @enum BoundaryConditionOption
   *
   * The options for essential boundary conditions
   */
  enum struct BoundaryConditionOption : integer
  {
    OUTFLOW,    //!<Outflow
    SYMMETRY    //!<Symmetry
  };

  /**
   * Constructor
   * @param name The name of the solver instance
   * @param parent the parent group of the solver
   */
  SolidMechanicsMPM( const string & name,
                     Group * const parent );


  SolidMechanicsMPM( SolidMechanicsMPM const & ) = delete;
  SolidMechanicsMPM( SolidMechanicsMPM && ) = default;

  SolidMechanicsMPM & operator=( SolidMechanicsMPM const & ) = delete;
  SolidMechanicsMPM & operator=( SolidMechanicsMPM && ) = delete;

  /**
   * destructor
   */
  virtual ~SolidMechanicsMPM() override;

  /**
   * @return The string that may be used to generate a new instance from the SolverBase::CatalogInterface::CatalogType
   */
  static string catalogName() { return "SolidMechanics_MPM"; }
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
  virtual
  real64 solverStep( real64 const & time_n,
                     real64 const & dt,
                     integer const cycleNumber,
                     DomainPartition & domain ) override;

  virtual
  real64 explicitStep( real64 const & time_n,
                       real64 const & dt,
                       integer const cycleNumber,
                       DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override final
  {
    // There should be nothing to update
    GEOS_UNUSED_VAR( domain );
  };

  /**@}*/

  /**
   * This method is called when its host event is triggered
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;


  template< typename CONSTITUTIVE_BASE,
            typename KERNEL_WRAPPER,
            typename ... PARAMS >
  void assemblyLaunch( DomainPartition & domain,
                       DofManager const & dofManager,
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs,
                       PARAMS && ... params );


  template< typename ... PARAMS >
  real64 explicitKernelDispatch( MeshLevel & mesh,
                                 arrayView1d< string const > const & targetRegions,
                                 string const & finiteElementName,
                                 real64 const dt,
                                 std::string const & elementListName );

  /**
   * Applies displacement boundary conditions to the system for implicit time integration
   * @param time The time to use for any lookups associated with this BC
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain The DomainPartition.
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   */

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr char const * cflFactorString() { return "cflFactor"; }
    static constexpr char const * timeIntegrationOptionString() { return "timeIntegrationOption"; }
    static constexpr char const * solidMaterialNamesString() { return "solidMaterialNames"; }
    static constexpr char const * forceExternalString() { return "externalForce"; }
    static constexpr char const * forceInternalString() { return "internalForce"; }
    static constexpr char const * massString() { return "mass"; }
    static constexpr char const * velocityString() { return "velocity"; }
    static constexpr char const * momentumString() { return "momentum"; }
    static constexpr char const * accelerationString() { return "acceleration"; }
    static constexpr char const * forceContactString() { return "contactForce"; }
    static constexpr char const * damageString() { return "damage"; }
    static constexpr char const * damageGradientString() { return "damageGradient"; }
    static constexpr char const * maxDamageString() { return "maxDamage"; }
    static constexpr char const * surfaceNormalString() { return "surfaceNormal"; }
    static constexpr char const * materialPositionString() { return "materialPosition"; }

    static constexpr char const * boundaryNodesString() { return "boundaryNodes"; }
    static constexpr char const * bufferNodesString() { return "bufferNodes"; }

    dataRepository::ViewKey timeIntegrationOption = { timeIntegrationOptionString() };
  } solidMechanicsViewKeys;

  void initialize( NodeManager & nodeManager,
                   ParticleManager & particleManager,
                   SpatialPartition & partition );

  void resizeGrid( SpatialPartition & partition,
                   NodeManager & nodeManager,
                   real64 const dt );

  void syncGridFields( std::vector< std::string > const & fieldNames,
                       DomainPartition & domain,
                       NodeManager & nodeManager,
                       MeshLevel & mesh,
                       MPI_Op op );

  void singleFaceVectorFieldSymmetryBC( const int face,
                                        arrayView3d< real64 > const & vectorMultiField,
                                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                        Group & nodeSets );

  void enforceGridVectorFieldSymmetryBC( arrayView3d< real64 > const & vectorMultiField,
                                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                         Group & nodeSets );

  void applyEssentialBCs( const real64 dt,
                          const real64 time_n,
                          NodeManager & nodeManager );

  void computeGridSurfaceNormals( ParticleManager & particleManager,
                                  NodeManager & nodeManager );

  void normalizeGridSurfaceNormals( arrayView2d< real64 const > const & gridMass,
                                    arrayView3d< real64 > const & gridSurfaceNormal );

  void computeContactForces( real64 const dt,
                             arrayView2d< real64 const > const & gridMass,
                             arrayView2d< real64 const > const & gridDamage,
                             arrayView2d< real64 const > const & gridMaxDamage,
                             arrayView3d< real64 const > const & gridVelocity,
                             arrayView3d< real64 const > const & gridMomentum,
                             arrayView3d< real64 const > const & gridSurfaceNormal,
                             arrayView3d< real64 const > const & gridMaterialPosition,
                             arrayView3d< real64 > const & gridContactForce );

  void computePairwiseNodalContactForce( int const & separable,
                                         real64 const & dt,
                                         real64 const & mA,
                                         real64 const & mB,
                                         arraySlice1d< real64 const > const vA,
                                         arraySlice1d< real64 const > const GEOS_UNUSED_PARAM( vB ),
                                         arraySlice1d< real64 const > const qA,
                                         arraySlice1d< real64 const > const qB,
                                         arraySlice1d< real64 const > const nA,
                                         arraySlice1d< real64 const > const nB,
                                         arraySlice1d< real64 const > const xA, // Position of field A
                                         arraySlice1d< real64 const > const xB, // Position of field B
                                         arraySlice1d< real64 > const fA,
                                         arraySlice1d< real64 > const fB );

  void computeOrthonormalBasis( const real64 * e1,  // input "normal" unit vector.
                                real64 * e2,        // output "tangential" unit vector.
                                real64 * e3 );      // output "tangential" unit vector.

  void setGridFieldLabels( NodeManager & nodeManager );

  void solverProfiling( std::string label );

  void solverProfilingIf( std::string label, bool condition );

  real64 computeNeighborList( ParticleManager & particleManager );

  void optimizeBinSort( ParticleManager & particleManager );

  real64 kernel( real64 const & r ); // distance from particle to query point

  void kernelGradient( arraySlice1d< real64 const > const x,  // query point
                       std::vector< real64 > & xp,            // particle location
                       real64 const & r,                      // distance from particle to query point
                       real64 * result );

  real64 computeKernelField( arraySlice1d< real64 const > const x,    // query point
                             arrayView2d< real64 const > const xp,    // List of neighbor particle locations.
                             arrayView1d< real64 const > const Vp,    // List of neighbor particle volumes.
                             arrayView1d< real64 const > const fp );  // scalar field values (e.g. damage) at neighbor particles

  void computeKernelFieldGradient( arraySlice1d< real64 const > const x,       // query point
                                   std::vector< std::vector< real64 > > & xp,  // List of neighbor particle locations.
                                   std::vector< real64 > & Vp,                 // List of neighbor particle volumes.
                                   std::vector< real64 > & fp,                 // scalar field values (e.g. damage) at neighbor particles
                                   arraySlice1d< real64 > const result );

  void computeKernelVectorGradient( arraySlice1d< real64 const > const x,       // query point
                                    std::vector< std::vector< real64 > > & xp,  // List of neighbor particle locations.
                                    std::vector< real64 > & Vp,                 // List of neighbor particle volumes.
                                    std::vector< std::vector< real64 > > & fp,  // vector field values (e.g. velocity) at neighbor particles
                                    arraySlice2d< real64 > const result );

  void computeDamageFieldGradient( ParticleManager & particleManager );

  void updateSurfaceFlagOverload( ParticleManager & particleManager );

  void projectDamageFieldGradientToGrid( ParticleManager & particleManager,
                                         NodeManager & nodeManager );

  void updateDeformationGradient( real64 dt,
                                  ParticleManager & particleManager );

  void updateConstitutiveModelDependencies( ParticleManager & particleManager );

  void updateStress( real64 dt,
                     ParticleManager & particleManager );

  void particleKinematicUpdate( ParticleManager & particleManager );

  void computeAndWriteBoxAverage( const real64 dt,
                                  const real64 time_n,
                                  ParticleManager & particleManager );

  void initializeGridFields( NodeManager & nodeManager );

  void boundaryConditionUpdate( real64 dt, real64 time_n );

  void particleToGrid( ParticleManager & particleManager,
                       NodeManager & nodeManager );

  void gridTrialUpdate( real64 dt,
                        NodeManager & nodeManager );

  void enforceContact( real64 dt,
                       DomainPartition & domain,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager,
                       MeshLevel & mesh );

  void interpolateFTable( real64 dt, real64 time_n );

  void gridToParticle( real64 dt,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager );

  void updateSolverDependencies( ParticleManager & particleManager );

  real64 getStableTimeStep( ParticleManager & particleManager );

  void deleteBadParticles( ParticleManager & particleManager );

  void printProfilingResults();

  void computeSurfaceFlags( ParticleManager & particleManager );

  void computeSphF( ParticleManager & particleManager );

  // void directionalOverlapCorrection( real64 dt, ParticleManager & particleManager );

  int evaluateSeparabilityCriterion( localIndex const & A,
                                     localIndex const & B,
                                     real64 const & damageA,
                                     real64 const & damageB,
                                     real64 const & maxDamageA,
                                     real64 const & maxDamageB );

  void flagOutOfRangeParticles( ParticleManager & particleManager );

  void computeRVectors( ParticleManager & particleManager );

  void cpdiDomainScaling( ParticleManager & particleManager );

  void resizeMappingArrays( ParticleManager & particleManager );

  void populateMappingArrays( ParticleManager & particleManager,
                              NodeManager & nodeManager );

protected:
  virtual void postInputInitialization() override final;

  std::vector< array2d< localIndex > > m_mappedNodes; // mappedNodes[subregion index][particle index][node index]. dims = {# of subregions,
                                                      // # of particles, # of nodes a particle on the subregion maps to}
  std::vector< array2d< real64 > > m_shapeFunctionValues; // mappedNodes[subregion][particle][nodal shape function value]. dims = {# of
                                                          // subregions, # of particles, # of nodes a particle on the subregion maps to}
  std::vector< array3d< real64 > > m_shapeFunctionGradientValues; // mappedNodes[subregion][particle][nodal shape function gradient
                                                                  // value][direction]. dims = {# of subregions, # of particles, # of nodes
                                                                  // a particle on the subregion maps to, 3}

  int m_solverProfiling;
  std::vector< real64 > m_profilingTimes;
  std::vector< std::string > m_profilingLabels;

  TimeIntegrationOption m_timeIntegrationOption;

  int m_prescribedBcTable;
  array1d< int > m_boundaryConditionTypes; // TODO: Surely there's a way to have just one variable here
  array2d< real64 > m_bcTable;

  int m_prescribedBoundaryFTable;
  Path m_fTablePath;
  int m_fTableInterpType;
  array2d< real64 > m_fTable;
  array1d< real64 > m_domainF;
  array1d< real64 > m_domainL;

  int m_boxAverageHistory;
  int m_reactionHistory;

  int m_needsNeighborList;
  real64 m_neighborRadius;
  int m_binSizeMultiplier;

  int m_useDamageAsSurfaceFlag;

  int m_cpdiDomainScaling;

  real64 m_smallMass;

  int m_numContactGroups, m_numContactFlags, m_numVelocityFields;
  real64 m_separabilityMinDamage;
  int m_treatFullyDamagedAsSingleField;
  int m_surfaceDetection;
  int m_damageFieldPartitioning;
  int m_contactGapCorrection;
  // int m_directionalOverlapCorrection;
  real64 m_frictionCoefficient;

  int m_planeStrain;
  int m_numDims;

  array1d< real64 > m_hEl;                // Grid spacing in x-y-z
  array1d< real64 > m_xLocalMin;          // Minimum local grid coordinate including ghost nodes
  array1d< real64 > m_xLocalMax;          // Maximum local grid coordinate including ghost nodes
  array1d< real64 > m_xLocalMinNoGhost;   // Minimum local grid coordinate EXCLUDING ghost nodes
  array1d< real64 > m_xLocalMaxNoGhost;   // Maximum local grid coordinate EXCLUDING ghost nodes
  array1d< real64 > m_xGlobalMin;         // Minimum global grid coordinate excluding buffer nodes
  array1d< real64 > m_xGlobalMax;         // Maximum global grid coordinate excluding buffer nodes
  array1d< real64 > m_partitionExtent;    // Length of each edge of partition including buffer and ghost cells
  array1d< real64 > m_domainExtent;       // Length of each edge of global domain excluding buffer cells
  array1d< int > m_nEl;                   // Number of elements in each grid direction including buffer and ghost cells
  array3d< int > m_ijkMap;        // Map from indices in each spatial dimension to local node ID

private:
  struct BinKey
  {
    localIndex regionIndex;
    localIndex subRegionIndex;
    localIndex binIndex;

    bool operator==( BinKey const & other ) const
    {
      return (regionIndex == other.regionIndex && subRegionIndex == other.subRegionIndex && binIndex == other.binIndex);
    }
  };

  struct BinKeyHash
  {
    std::size_t operator()( BinKey const & k ) const
    {
      using std::size_t;
      using std::hash;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:
      return ((std::hash< localIndex >()( k.regionIndex )
               ^ (std::hash< localIndex >()( k.subRegionIndex ) << 1)) >> 1)
             ^ (std::hash< localIndex >()( k.binIndex ) << 1);
    }
  };

  void setParticlesConstitutiveNames( ParticleSubRegionBase & subRegion ) const;
};

ENUM_STRINGS( SolidMechanicsMPM::TimeIntegrationOption,
              "QuasiStatic",
              "ImplicitDynamic",
              "ExplicitDynamic" );

ENUM_STRINGS( SolidMechanicsMPM::BoundaryConditionOption,
              "Outflow",
              "Symmetry" );

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************


} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_ */
