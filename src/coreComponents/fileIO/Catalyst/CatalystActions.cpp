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
 * @file CatalystActions.cpp
 */

// source includes
#include "CatalystActions.hpp"

// external includes
#include <catalyst.hpp>

#include <functional>

namespace geos
{

namespace internal
{

/*
 * Convert the encapsulated simulation conduit node into a catalyst conduit node
 */
conduit_cpp::Node Convert(ConduitCapsule* simulationNode)
{
  conduit_cpp::Node convertedNode;

  if (!simulationNode)
  {
    return convertedNode;
  }

  if(simulationNode->IsList() && simulationNode->IsInt8())
  {
    convertedNode.set_int8_vector(simulationNode->GetInt8Array());
  }
  else if(simulationNode->IsList() && simulationNode->IsInt16())
  {
    convertedNode.set_int16_vector(simulationNode->GetInt16Array());
  }
  else if(simulationNode->IsList() && simulationNode->IsInt32())
  {
      convertedNode.set_int32_vector(simulationNode->GetInt32Array());
  }
  else if(simulationNode->IsList() && simulationNode->IsInt64())
  {
      convertedNode.set_int64_vector(simulationNode->GetInt64Array());
  }
  else if(simulationNode->IsList() && simulationNode->IsUInt8())
  {
    convertedNode.set_uint8_vector(simulationNode->GetUInt8Array());
  }
  else if(simulationNode->IsList() && simulationNode->IsUInt16())
  {
    convertedNode.set_uint16_vector(simulationNode->GetUInt16Array());
  }
  else if(simulationNode->IsList() && simulationNode->IsUInt32())
  {
      convertedNode.set_uint32_vector(simulationNode->GetUInt32Array());
  }
  else if(simulationNode->IsList() && simulationNode->IsUInt64())
  {
      convertedNode.set_uint64_vector(simulationNode->GetUInt64Array());
  }
  else if(simulationNode->IsList() && simulationNode->IsFloat32())
  {
      convertedNode.set_float32_vector(simulationNode->GetFloat32Array());
  }
  else if(simulationNode->IsList() && simulationNode->IsFloat64())
  {
      convertedNode.set_float64_vector(simulationNode->GetFloat64Array());
  }
  else if(simulationNode->IsInt8())
  {
    convertedNode.set_int8(simulationNode->GetInt8Value());
  }
  else if(simulationNode->IsInt16())
  {
    convertedNode.set_int16(simulationNode->GetInt16Value());
  }
  else if(simulationNode->IsInt32())
  {
      convertedNode.set_int32(simulationNode->GetInt32Value());
  }
  else if(simulationNode->IsInt64())
  {
      convertedNode.set_int64(simulationNode->GetInt64Value());
  }
  else if(simulationNode->IsUInt8())
  {
    convertedNode.set_uint8(simulationNode->GetUInt8Value());
  }
  else if(simulationNode->IsUInt16())
  {
    convertedNode.set_uint16(simulationNode->GetUInt16Value());
  }
  else if(simulationNode->IsUInt32())
  {
      convertedNode.set_uint32(simulationNode->GetUInt32Value());
  }
  else if(simulationNode->IsUInt64())
  {
      convertedNode.set_uint64(simulationNode->GetUInt64Value());
  }
  else if(simulationNode->IsFloat32())
  {
      convertedNode.set_float32(simulationNode->GetFloat32Value());
  }
  else if(simulationNode->IsFloat64())
  {
      convertedNode.set_float64(simulationNode->GetFloat64Value());
  }
  else if(simulationNode->IsString())
  {
    convertedNode.set_string(simulationNode->GetStringValue());
  }

  // Recursively convert children
  for(std::size_t iChild = 0; iChild < simulationNode->GetNumberOfChildren(); ++iChild)
  {
    auto childConduitNode = simulationNode->GetChild(iChild);
    auto name = childConduitNode->GetName();
    if(name.empty())
    {
      name="root";
    }
    convertedNode.set_path_node(name, Convert(childConduitNode.get()));
  }

  return convertedNode;
}

bool CatalystAction(ConduitCapsule* simulationNode, std::function<bool(conduit_node_impl*)> action)
{
  if(!simulationNode)
  {
    return false;
  }

  conduit_cpp::Node convertedNode = Convert(simulationNode);

  return (action(conduit_cpp::c_node(&convertedNode)) == catalyst_status::catalyst_status_ok);
}

}

bool CatalystInitialize(ConduitCapsule* simulationNode)
{
  return internal::CatalystAction(simulationNode, catalyst_initialize);
}

bool CatalystExecute(ConduitCapsule* simulationNode)
{
  return internal::CatalystAction(simulationNode, catalyst_execute);
}

bool CatalystFinalize(ConduitCapsule* simulationNode)
{
  return internal::CatalystAction(simulationNode, catalyst_finalize);
}

}
