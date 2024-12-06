/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior
 * University Copyright (c) 2018-2020 TotalEnergies Copyright (c) 2019- GEOSX
 * Contributors All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS
 * files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CatalystOutput.cpp
 */

// GEOS
#include "CatalystOutput.hpp"

#if defined(GEOSX_USE_PYGEOSX)
#include "fileIO/python/PyCatalystOutputType.hpp"
#endif

#include "fileIO/Catalyst/CatalystActions.hpp"
#include "fileIO/Catalyst/GenericConduitCapsule.tpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

// TPL
#include <conduit.hpp>
#include <conduit_cpp_to_c.hpp>

// std
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void SanitizeSimpleNode(conduit::Node& channel) {
  auto& blueprintNode = channel["data"];
  if (blueprintNode["topologies"].number_of_children() > 1) {
    GEOS_ERROR(
        "CatalystOutput: multi region passed to simple region code path");
    return;
  }
  std::string topoName = blueprintNode["topologies"].child(0).name();
  auto& fields = blueprintNode["fields"];
  for (conduit::index_t iField = 0; iField < fields.number_of_children();
       ++iField) {
    auto& field = fields.child(iField);
    if (field["topology"].as_char8_str() == std::string("nodes")) {
      field["topology"] = topoName;
      field["association"] = "vertex";
    }
  }
}

void SanitizeMultiNodeInSitu(conduit::Node& channel) {
  auto& blueprintNode = channel["data"];
  for (conduit::index_t iTopo = 0; iTopo < blueprintNode["topologies"].number_of_children(); ++iTopo) {
    auto& topology = blueprintNode["topologies"].child(iTopo);
    std::string topoName = topology.name();
    auto& iMesh = blueprintNode[topoName];

    iMesh["state"].update(blueprintNode["state"]);

    iMesh["coordsets"].update(blueprintNode["coordsets"]);
    iMesh["topologies/" + topoName].update(topology);
    iMesh["fields"].update(topology["fields"]);

    // Deprec ?
    auto& dstFields = iMesh["fields"];
    const auto& srcFields = blueprintNode["fields"];
    for (conduit::index_t iField = 0; iField < srcFields.number_of_children(); ++iField) {
      const auto& srcField = srcFields.child(iField);
      if (srcField["topology"].as_char8_str() == std::string("nodes")) {
        auto& dstField = dstFields[srcField.name()];
        dstField.update(srcField);
        dstField["topology"] = topoName;
        dstField["association"] = "vertex";
      } else if (srcField["topology"].as_char8_str() == topoName) {
        auto& dstField = dstFields[srcField.name()];
        dstField.update(srcField);
      }
    }
    if (dstFields.number_of_children() == 0) {
      iMesh.remove("fields");
    }
  }
  blueprintNode.remove("state");
  blueprintNode.remove("topologies");
  blueprintNode.remove("coordsets");
  blueprintNode.remove("fields");
}

struct ScopedData {};

template <typename T>
struct ScopedTypedData : public ScopedData {
 public:
  ScopedTypedData(T data) : Data(data) {}

 private:
  T Data;
};

struct ScopedDataContainer : public std::vector<std::unique_ptr<ScopedData>> {
 public:
  template <typename T>
  void Add(T data) {
    this->emplace_back(
        std::unique_ptr<ScopedData>(new ScopedTypedData<T>(data)));
  }
};

template <typename IndexT>
std::set<std::uint64_t> ConstructNodeMapping(const IndexT* connectivity,
                                             conduit::index_t size) {
  std::set<std::uint64_t> nodeMapping;
  for (conduit::index_t iConnect = 0; iConnect < size; ++iConnect) {
    nodeMapping.emplace(
        static_cast<const std::uint64_t>(connectivity[iConnect]));
  }
  return nodeMapping;
}

std::set<std::uint64_t> ConstructNodeMapping(
    const conduit::Node& connectivity) {
  std::string dtype = connectivity.dtype().name();
  conduit::index_t numValues = connectivity.dtype().number_of_elements();
#define DispatchNameToType(name, type)                                   \
  do {                                                                   \
    if (dtype == name) {                                                 \
      auto res = ConstructNodeMapping(                                   \
          static_cast<const type*>(connectivity.data_ptr()), numValues); \
      return res;                                                        \
    }                                                                    \
  } while (0)
  DispatchNameToType("int8", signed char);
  DispatchNameToType("int16", short);
  DispatchNameToType("int32", std::int32_t);
  DispatchNameToType("int64", std::int64_t);
  DispatchNameToType("uint8", unsigned char);
  DispatchNameToType("uint16", unsigned short);
  DispatchNameToType("uint32", std::uint32_t);
  DispatchNameToType("uint64", std::uint64_t);
  DispatchNameToType("float32", float);
  DispatchNameToType("float64", double);
#undef DispatchNameToType
  GEOS_ERROR("Could not dispatch connectivity with type " + dtype);
  return std::set<std::uint64_t>();
}

template <typename DataT>
std::shared_ptr<std::vector<DataT>> ConstructConnectivity(
    const std::set<std::uint64_t>& nodeMapping, const DataT* src,
    conduit_index_t size) {
  auto res = std::make_shared<std::vector<DataT>>(size);
  for (conduit::index_t iConnect = 0; iConnect < size; ++iConnect) {
    (*res)[iConnect] = std::distance(
        nodeMapping.begin(),
        nodeMapping.find(static_cast<std::uint64_t>(src[iConnect])));
  }
  return res;
}

bool DispatchConnectivity(const std::set<std::uint64_t>& nodeMapping,
                          const conduit::Node& src, conduit::Node& dst,
                          ScopedDataContainer& dataContainer) {
  std::string dtype = src.dtype().name();
  conduit_index_t size = src.dtype().number_of_elements();
#define DispatchNameToType(name, type)                                  \
  do {                                                                  \
    if (dtype == name) {                                                \
      auto res = ConstructConnectivity(                                 \
          nodeMapping, static_cast<const type*>(src.data_ptr()), size); \
      dataContainer.Add(res);                                           \
      dst.set_external(res->data(), res->size());                       \
      return true;                                                      \
    }                                                                   \
  } while (0)
  DispatchNameToType("int8", signed char);
  DispatchNameToType("int16", short);
  DispatchNameToType("int32", std::int32_t);
  DispatchNameToType("int64", std::int64_t);
  DispatchNameToType("uint8", unsigned char);
  DispatchNameToType("uint16", unsigned short);
  DispatchNameToType("uint32", std::uint32_t);
  DispatchNameToType("uint64", std::uint64_t);
  DispatchNameToType("float32", float);
  DispatchNameToType("float64", double);
#undef DispatchNameToType
  GEOS_ERROR("Could not dispatch connectivity type " + dtype);
  return false;
}

template <typename DataT>
std::shared_ptr<std::vector<DataT>> ConstructReducedField(
    const std::set<std::uint64_t>& nodeMapping, const DataT* src,
    std::size_t offset, std::size_t stride) {
  auto res = std::make_shared<std::vector<DataT>>(nodeMapping.size(), 0);
  auto resIt = res->begin();
  for (const auto& nodeIndex : nodeMapping) {
    *resIt = src[static_cast<std::size_t>(nodeIndex) * stride + offset];
    resIt++;
  }
  return res;
}

bool DispatchFieldToContainers(const std::set<std::uint64_t>& nodeMapping,
                               const conduit::Node& src, conduit::Node& dst,
                               ScopedDataContainer& dataContainer,
                               std::size_t offset = 0, std::size_t stride = 0) {
  std::string dtype = src.dtype().name();
#define DispatchNameToType(name, type)                                   \
  do {                                                                   \
    if (dtype == name) {                                                 \
      auto res = ConstructReducedField(                                  \
          nodeMapping, static_cast<const type*>(src.data_ptr()), offset, \
          stride);                                                       \
      dataContainer.Add(res);                                            \
      dst.set_external(res->data(), res->size());                        \
      return true;                                                       \
    }                                                                    \
  } while (0)
  DispatchNameToType("int8", signed char);
  DispatchNameToType("int16", short);
  DispatchNameToType("int32", std::int32_t);
  DispatchNameToType("int64", std::int64_t);
  DispatchNameToType("uint8", unsigned char);
  DispatchNameToType("uint16", unsigned short);
  DispatchNameToType("uint32", std::uint32_t);
  DispatchNameToType("uint64", std::uint64_t);
  DispatchNameToType("float32", float);
  DispatchNameToType("float64", double);
#undef DispatchNameToType
  GEOS_ERROR("Could not dispatch field type " + dtype);
  return false;
}

void SanitizeMultiNodeInTransit(conduit::Node& channel,
                                ScopedDataContainer& dataContainer) {
  auto& blueprintNode = channel["data"];
  for (conduit::index_t iTopo = 0;
       iTopo < blueprintNode["topologies"].number_of_children(); ++iTopo) {
    const auto& topology = blueprintNode["topologies"].child(iTopo);

    std::string topoName = topology.name();
    auto& iMesh = blueprintNode[topoName];

    iMesh["state"].update(blueprintNode["state"]);

    const auto& srcConnectivity = topology["elements/connectivity"];
    std::set<std::uint64_t> nodeMapping = ConstructNodeMapping(srcConnectivity);

    iMesh["topologies/" + topoName + "/type"].update(topology["type"]);
    iMesh["topologies/" + topoName + "/elements/shape"].update(
        topology["elements/shape"]);
    iMesh["topologies/" + topoName + "/coordset"].update(topology["coordset"]);

    auto& dstConnectivity =
        iMesh["topologies/" + topoName + "/elements/connectivity"];
    if (!DispatchConnectivity(nodeMapping, srcConnectivity, dstConnectivity,
                              dataContainer)) {
      GEOS_ERROR(
          "CatalystOutput: SanitizeMultiNodeInTransit: Dispatch of "
          "connecitivities failed");
    }

    std::string coordsetName =
        std::string("coordsets/") + topology["coordset"].as_char8_str();
    auto& coords = blueprintNode[coordsetName];
    iMesh[coordsetName + "/type"].update(coords["type"]);
    std::map<std::string, const conduit::Node&> xyz = {
        {"x", coords["values/x"]},
        {"y", coords["values/y"]},
        {"z", coords["values/z"]}};
    std::size_t offset = 0;
    std::size_t stride = 1;
    if (!coords.is_contiguous()) {
      stride = 3;
    }
    for (const auto& component : xyz) {
      auto& dst = iMesh[coordsetName + "/values/" + component.first];
      if (!DispatchFieldToContainers(nodeMapping, component.second, dst,
                                     dataContainer, offset, stride)) {
        GEOS_ERROR(
            "CatalystOutput: SanitizeMultiNodeInTransit: Could not "
            "dispatch " +
            component.second.name() + " field.");
      }
    }

    const auto& srcFields = blueprintNode["fields"];
    auto& dstFields = iMesh["fields"];
    for (conduit::index_t iField = 0; iField < srcFields.number_of_children();
         ++iField) {
      const auto& srcField = srcFields.child(iField);
      bool isPointData =
          srcField["topology"].as_char8_str() == std::string("nodes");
      bool isCellData = (srcField["topology"].as_char8_str() == topoName);
      if (!isPointData && !isCellData) {
        continue;
      }
      auto& dstField = dstFields[srcField.name()];
      if (isCellData) {
        dstField.update(srcField);
        continue;
      }
      dstField["topology"] = topoName;
      dstField["volume_dependent"].update(srcField["volume_dependent"]);
      dstField["association"] = "vertex";
      if (!DispatchFieldToContainers(nodeMapping, srcField["values"],
                                     dstField["values"], dataContainer)) {
        GEOS_ERROR(
            "CatalystOutput: SanitizeMultiNodeInTransit: Could not "
            "dispatch " +
            srcField.name() + " field.");
      }
    }
    if (dstFields.number_of_children() == 0) {
      iMesh.remove("fields");
    }
  }
  blueprintNode.remove("state");
  blueprintNode.remove("topologies");
  blueprintNode.remove("coordsets");
  blueprintNode.remove("fields");
}

void SanitizeNode(conduit::Node& channel, bool isInTransit,
                  ScopedDataContainer& dataContainer) {
  channel.remove("data/topologies/nodes");
  if (channel["data/topologies"].number_of_children() > 1) {
    channel["type"] = "multimesh";
    if (isInTransit) {
      SanitizeMultiNodeInTransit(channel, dataContainer);
      return;
    }
    channel.remove("data/fields"); // fields info to use is inside each topologies
    SanitizeMultiNodeInSitu(channel);
    return;
  }

  SanitizeSimpleNode(channel);
  return;
}

}  // namespace

namespace geos {

struct CatalystOutput::CatalystInternals {
  bool initialized = false;
  std::string scripts = "";
  std::string implementation = "";
  std::string implementationPath = "";
  std::string adiosConfig = "";
  std::string channelName = "fullfield";
  std::string sstFileName = "gs.bp";

  void initializeCatalyst() {
    std::string scriptsCopy = this->scripts;
    std::vector<std::string> scriptList;
    std::string delimiter = ":";
    std::size_t pos = scriptsCopy.find(delimiter);
    scriptList.emplace_back(scriptsCopy.substr(0, pos));
    scriptsCopy.erase(0, pos + delimiter.length());
    while ((pos = scriptsCopy.find(delimiter)) != std::string::npos) {
      scriptList.emplace_back(scriptsCopy.substr(0, pos));
      scriptsCopy.erase(0, pos + delimiter.length());
    }

    if (scriptList.empty() || scriptList[0] == "") {
      GEOS_ERROR("CatalystOutput: Constructor: no catalyst scripts found.");
    }

    conduit::Node initializer;
    for (std::size_t iScr = 0; iScr < scriptList.size(); ++iScr) {
      initializer["catalyst/scripts/script" + std::to_string(iScr)] =
          scriptList[iScr];
    }

    if (!this->implementation.empty()) {
      initializer["catalyst_load/implementation"] = this->implementation;
      if (!this->implementationPath.empty()) {
        initializer["catalyst_load/search_paths/" + this->implementation] =
            this->implementationPath;
      }
    }

    if (this->implementation == "adios" ||
        (this->implementation.empty() &&
         getEnv("CATALYST_IMPLEMENTATION_NAME") == "adios")) {
      if (this->adiosConfig == "") {
        GEOS_ERROR(
            "CatalystOutput: Constructor: Cannot use the catalyst/adios2 "
            "implementation without providing an adios2 configuration file.");
      }
      initializer["adios/config"] = this->adiosConfig;
    }

    std::string envSSTfilename = getEnv("CATALYST_SST_FILTENAME");
    if( envSSTfilename != "") {
      this->sstFileName = envSSTfilename;
    }

    auto capsule = GenericConduitCapsule<conduit::Node>(&initializer);
    if (!CatalystInitialize(&capsule)) {
      GEOS_ERROR("CatalystOutput: Constructor: catalyst failed to initialize.");
    }

    this->initialized = true;
  }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
CatalystOutput::CatalystOutput(std::string const& name, Group* const parent)
    : BlueprintOutput(name, parent),
      internal(std::unique_ptr<CatalystInternals>(new CatalystInternals())) {
  this->registerWrapper("scripts", &this->internal->scripts)
      .setApplyDefaultValue("")
      .setInputFlag(dataRepository::InputFlags::REQUIRED)
      .setDescription("Column separated paths to the catalyst scripts.");

  this->registerWrapper("implementation", &this->internal->implementation)
      .setApplyDefaultValue("")
      .setInputFlag(dataRepository::InputFlags::OPTIONAL)
      .setDescription("Name of the catalyst implementation to use.");

  this->registerWrapper("implementationPath",
                        &this->internal->implementationPath)
      .setApplyDefaultValue("")
      .setInputFlag(dataRepository::InputFlags::OPTIONAL)
      .setDescription("Path to the catalyst the implementation to use.");

  this->registerWrapper("adiosConfig", &this->internal->adiosConfig)
      .setApplyDefaultValue("")
      .setInputFlag(dataRepository::InputFlags::OPTIONAL)
      .setDescription(
          "Path to the adios configuration file when using the catalyst-adios "
          "implementation.");

  this->registerWrapper("fullFieldChannelName", &this->internal->channelName)
      .setApplyDefaultValue("")
      .setInputFlag(dataRepository::InputFlags::OPTIONAL)
      .setDescription(
          "Name to give to the channel passing the full field data.");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
CatalystOutput::~CatalystOutput() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
bool CatalystOutput::execute(real64 const time_n, real64 const /*dt*/,
                             integer const cycleNumber,
                             integer const /*eventCounter*/,
                             real64 const /*eventProgress*/,
                             DomainPartition& domain) {
  GEOS_MARK_FUNCTION;

  if (!this->internal->initialized) {
    this->internal->initializeCatalyst();
  }

  conduit::Node executeRoot;
  auto& catalystState = executeRoot["catalyst/state"];
  catalystState["timestep"].set(cycleNumber);
  catalystState["time"].set(time_n);
  catalystState["sstFileName"].set(this->internal->sstFileName);

  auto& channel =
      executeRoot["catalyst/channels/" + this->internal->channelName];
  channel["type"] = "mesh";

  auto& meshGEOSRoot = channel["data"];
  this->mapMesh(time_n, cycleNumber, domain, meshGEOSRoot);

  bool isInTransit =
      this->internal->implementation == "adios" ||
      (this->internal->implementation.empty() &&
       getEnv("CATALYST_IMPLEMENTATION_NAME") == "adios");

  ScopedDataContainer dataScoping;
  ::SanitizeNode(channel, isInTransit, dataScoping);

  auto capsule = GenericConduitCapsule<conduit::Node>(&executeRoot);
  if (!CatalystExecute(&capsule)) {
    GEOS_ERROR("CatalystOutput: execute: catalyst failed to execute.");
    return true;
  }

  return false;
}

const std::string CatalystOutput::getEnv( const char* varname )
{
  char* varptr = std::getenv( varname );
  return varptr == nullptr ? "" : std::string( varptr );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void CatalystOutput::cleanup(real64 const time_n, integer const cycleNumber,
                             integer const eventCounter,
                             real64 const eventProgress,
                             DomainPartition& domain) {
  GEOS_MARK_FUNCTION;

  if (this->execute(time_n, 0, cycleNumber, eventCounter, eventProgress,
                    domain)) {
    GEOS_ERROR("CatalystOutput: cleanup: last execute failed.");
  }

  conduit::Node emptyNode;
  auto capsule = GenericConduitCapsule<conduit::Node>(&emptyNode);
  if (!CatalystFinalize(&capsule)) {
    GEOS_ERROR("CatalystOutput: cleanup: finalize catalyst failed");
  }

  this->internal->initialized = false;
}

#if defined(GEOSX_USE_PYGEOSX)
PyTypeObject * CatalystOutput::getPythonType() const
{
  return python::getPyCatalystOutputType();
}
#endif

REGISTER_CATALOG_ENTRY(OutputBase, CatalystOutput, string const&,
                       dataRepository::Group* const)

}  // namespace geos
