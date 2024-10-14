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
 * @file BlueprintOutput.hpp
 */

#ifndef GEOS_FILEIO_OUTPUTS_BLUEPRINTOUTPUT_HPP_
#define GEOS_FILEIO_OUTPUTS_BLUEPRINTOUTPUT_HPP_

#include "fileIO/Outputs/OutputBase.hpp"

namespace geos
{

// Forward declarations
class MeshLevel;
class NodeManager;
class ElementRegionManager;

/**
 * @class BlueprintOutput
 * @brief A class for creating Conduit blueprint-based outputs.
 */
class BlueprintOutput : public OutputBase
{
public:

  /**
   * @brief Construct a new BlueprintOutput object.
   * @param name The name of the BlueprintObject in the data repository.
   * @param parent The parent Group.
   */
  BlueprintOutput( string const & name,
                   Group * const parent );

  /**
   * @brief Destructor.
   */
  virtual ~BlueprintOutput() override
  {}

  /**
   * @brief Get the name used to register this object in an XML file.
   * @return The string "Blueprint".
   */
  static string catalogName() { return "Blueprint"; }

  /**
   * @brief Writes out a Blueprint plot file.
   * @copydetails EventBase::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @brief Writes out a Blueprint plot file at the end of the simulation.
   * @copydetails ExecutableGroup::cleanup()
   */
  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override
  { execute( time_n, 0, cycleNumber, eventCounter, eventProgress, domain ); }

protected:

  /**
   * @brief Map the mesh to the conduit blueprint.
   * @param time The time value associated with the mesh
   * @param cycle The cycle from the execute method.
   * @param domain the domain from the execute method
   * @param meshRoot the conduit node to map the mesh to.
   */
  virtual void mapMesh( real64 const time,
                integer const cycle,
                DomainPartition & domain,
                conduit::Node & meshRoot );

private:

  /**
   * @brief Create the Blueprint coordinate set, the nodal topology and register nodal fields.
   * @param nodeManager The NodeManager to write out.
   * @param coordset The Blueprint coordinate set to populate.
   * @param topologies The Node containing all the topologies.
   * @param fields The Node containing all the fields.
   */
  void addNodalData( NodeManager const & nodeManager,
                     conduit::Node & coordset,
                     conduit::Node & topologies,
                     conduit::Node & fields );

  /**
   * @brief Create the blueprint topologies and register element fields.
   * @param elemRegionManager The ElementRegionManager to write out.
   * @param coordset The Blueprint coordinate set to associate the topologies with.
   * @param topologies The Node containing all the topologies.
   * @param fields The Node containing all the fields.
   */
  void addElementData( ElementRegionManager const & elemRegionManager,
                       conduit::Node & coordset,
                       conduit::Node & topologies,
                       conduit::Node & fields,
                       dataRepository::Group & averagedElementData );

  /**
   * @brief Write out all the children with the appropriate plot level of @p group.
   * @param group The group that contains the children to write out.
   * @param fields The Blueprint "fields" Node.
   * @param topology The name of the Blueprint "topology" that @p group is associated with.
   * @param prefix The string to prepend to the name of the field. If not specified no prefix is used.
   */
  void writeOutWrappersAsFields( Group const & group,
                                 conduit::Node & fields,
                                 string const & topology,
                                 string const & prefix="" );

  void writeOutConstitutiveData( dataRepository::Group const & constitutiveModel,
                                 conduit::Node & fields,
                                 string const & topology,
                                 dataRepository::Group & averagedElementData );

  // Used to determine which fields to write out.
  dataRepository::PlotLevel m_plotLevel = dataRepository::PlotLevel::LEVEL_1;

  // If true will write out the full quadrature data, otherwise it is averaged over.
  int m_outputFullQuadratureData = 0;
};


} // namespace geos

#endif // GEOS_FILEIO_OUTPUTS_BLUEPRINTOUTPUT_HPP_
