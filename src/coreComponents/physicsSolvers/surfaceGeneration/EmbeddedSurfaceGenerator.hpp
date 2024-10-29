/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EmbeddedSurfaceGenerator.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SURFACEGENERATION_EMBEDDEDSURFACEGENERATOR_HPP_
#define GEOS_PHYSICSSOLVERS_SURFACEGENERATION_EMBEDDEDSURFACEGENERATOR_HPP_

#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "mesh/DomainPartition.hpp"


namespace geos
{

struct NewObjectLists
{
  std::set< localIndex > newNodes;
  std::set< localIndex > newEdges;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > newElements;

  void insert( NewObjectLists const & lists );
};


class SpatialPartition;

class NodeManager;
class FaceManager;
class ElementRegionManager;
class ElementRegionBase;

/**
 * @class SurfaceGenerator
 *
 * This solver manages the mesh topology splitting methods.
 *
 */
class EmbeddedSurfaceGenerator : public SolverBase
{
public:
  EmbeddedSurfaceGenerator( const string & name,
                            Group * const parent );
  ~EmbeddedSurfaceGenerator() override;


  static string catalogName() { return "EmbeddedSurfaceGenerator"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( Group & MeshBody ) override final;

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const GEOS_UNUSED_PARAM( eventCounter ),
                        real64 const GEOS_UNUSED_PARAM( eventProgress ),
                        DomainPartition & domain ) override
  {
    solverStep( time_n, dt, cycleNumber, domain );
    return false;
  }

  /**
   * @brief xxx
   * @param[in] ...
   * @param[in] ...
   * @param[in] ...
   * @return ...
   *
   */
  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  /**@}*/

protected:

  virtual void initializePostSubGroups() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

  virtual void postRestartInitialization() override final
  {
    GEOS_ERROR( "Restarting is not supported for cases involving EmbeddedSurfaceGenerator" );
  }

private:

  void addToFractureStencil( DomainPartition & domain );

  void setGlobalIndices( ElementRegionManager & elemManager,
                         EmbeddedSurfaceNodeManager & embSurfNodeManager,
                         EmbeddedSurfaceSubRegion & embeddedSurfaceSubregion );

  void addEmbeddedElementsToSets( ElementRegionManager const & elemManager,
                                  EmbeddedSurfaceSubRegion & embeddedSurfaceSubregion );

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * solidMaterialNameString() {return "solidMaterialNames"; }
    constexpr static char const * fractureRegionNameString() {return "fractureRegion"; }
    constexpr static char const * targetObjectsNameString() {return "targetObjects"; }
    constexpr static char const * mpiCommOrderString() { return "mpiCommOrder"; }

    //TODO: rock toughness should be a material parameter, and we need to make rock toughness to KIC a constitutive
    // relation.
    constexpr static char const * rockToughnessString() {return "rockToughness"; }
  };

  // fracture region name
  string m_fractureRegionName;
  // target geometric objects to turn into fractures
  array1d< string > m_targetObjectsName;
  // Flag for consistent communication ordering
  int m_mpiCommOrder;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SURFACEGENERATION_EMBEDDEDSURFACEGENERATOR_HPP_ */
