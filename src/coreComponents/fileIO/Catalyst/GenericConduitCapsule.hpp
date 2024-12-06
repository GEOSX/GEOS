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
 * @file GenericConduitCapsule.hpp
 */

#ifndef GEOS_FILEIO_CATALYST_GENERICCONDUITCAPSULE_HPP_
#define GEOS_FILEIO_CATALYST_GENERICCONDUITCAPSULE_HPP_

// source includes
#include "ConduitCapsule.hpp"

namespace geos
{

/**
 * @class GenericConduitCapsule
 * @brief Helper template class for streamlining capsule creation for classic conduit nodes
 */
template<typename ConduitNodeT>
class GenericConduitCapsule : public ConduitCapsule
{
public:

  /**
   * @brief Constructor following RAII principles
   */
  GenericConduitCapsule(ConduitNodeT* node) : Node(node) {}

  /**
   * @brief Default destructor since there is no explicit memory responsibility
   */
  virtual ~GenericConduitCapsule() = default;

  /**
   * @brief Get the number of children this node element has in the hiearchy
   * @returns the number of children nodes this node has
   */
  virtual std::size_t GetNumberOfChildren() const;

  /**
   * @brief Get the iChild'th node below this node
   * @param iChild the index of the child node to get
   * @returns A shared pointer to the iChild node already encapsulated
   */
  virtual std::shared_ptr<ConduitCapsule> GetChild(std::size_t iChild);

  /**
   * @brief Get the name of the underlying node
   * @returns A string with the name of the node
   */
  virtual std::string GetName() const;

  ///@{
  /**
   * @brief Type checkers for the data held by the node
   * @returns true if the data type corresponds to the call
   */
  virtual bool IsList() const;
  virtual bool IsInt8() const;
  virtual bool IsInt16() const;
  virtual bool IsInt32() const;
  virtual bool IsInt64() const;
  virtual bool IsUInt8() const;
  virtual bool IsUInt16() const;
  virtual bool IsUInt32() const;
  virtual bool IsUInt64() const;
  virtual bool IsFloat32() const;
  virtual bool IsFloat64() const;
  virtual bool IsString() const;
  ///@}

  ///@{
  /**
   * @brief Check if the data held by the node is external
   */
  virtual bool IsDataExternal() const;
  /**
   * @brief Get a void* pointer to the external data when applicable
   * @returns pointer to external data when relevant and nullptr if not
   */
  virtual void* GetExternalDataPointer() const;
  /**
   * @brief Get the offset with which to read the external data
   * @retuns the offset meta-data
   */
  virtual long GetExternalDataOffset() const;
  /**
   * @brief Get the stride with which to read the external data
   * @retuns the stride meta-data
   */
  virtual long GetExternalDataStride() const;
  /**
   * @brief Get the number of bytes one element of external data occupies
   * @retuns the element size in bytes meta-data
   */
  virtual long GetExternalDataElementBytes() const;
  /**
   * @brief Get the endianness with which to read the external data
   * @retuns the endianness meta-data
   */
  virtual long GetExternalDataEndianness() const;
  ///@}


  ///@{
  /**
   * Get the values of data for conventional data types
   */
  virtual int8_t GetInt8Value() const;
  virtual int16_t GetInt16Value()  const;
  virtual int32_t GetInt32Value()  const;
  virtual int64_t GetInt64Value()  const;
  virtual uint8_t GetUInt8Value()   const;
  virtual uint16_t GetUInt16Value()  const;
  virtual uint32_t GetUInt32Value()  const;
  virtual uint64_t GetUInt64Value()  const;
  virtual float GetFloat32Value() const;
  virtual double GetFloat64Value() const;
  virtual std::string GetStringValue() const;
  ///@}

  ///@{
  /**
   * Get the value arrays of data for conventional data types
   */
  virtual std::vector<int8_t> GetInt8Array() const;
  virtual std::vector<int16_t> GetInt16Array() const;
  virtual std::vector<int32_t> GetInt32Array() const;
  virtual std::vector<int64_t> GetInt64Array() const;
  virtual std::vector<uint8_t> GetUInt8Array()  const;
  virtual std::vector<uint16_t> GetUInt16Array() const;
  virtual std::vector<uint32_t> GetUInt32Array() const;
  virtual std::vector<uint64_t> GetUInt64Array() const;
  virtual std::vector<float> GetFloat32Array() const;
  virtual std::vector<double> GetFloat64Array() const;
  ///@}

protected:

  /**
   * @brief the pointer to the node being encapsulated
   */
  ConduitNodeT* Node;

};

}

#endif
