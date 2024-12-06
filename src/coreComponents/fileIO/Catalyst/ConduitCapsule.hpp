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
 * @file ConduitCapsule.hpp
 */

#ifndef GEOS_FILEIO_CATALYST_CONDUITCAPSULE_HPP_
#define GEOS_FILEIO_CATALYST_CONDUITCAPSULE_HPP_

// External includes
#include <memory>
#include <string>
#include <vector>

namespace geos
{

/**
 * @class ConduitCapsule
 * @brief Pure interface for encapsulating a conduit node with the goal of coverting it to a catalyst conduit node
 */
class ConduitCapsule
{
public:

  /**
   * @brief No explicit memory deallocation here
   */
  virtual ~ConduitCapsule() = default;

  /**
   * @brief Get the number of children this node element has in the hiearchy
   * @returns the number of children nodes this node has
   */
  virtual std::size_t GetNumberOfChildren() const = 0;

  /**
   * @brief Get the iChild'th node below this node
   * @param iChild the index of the child node to get
   * @returns A shared pointer to the iChild node already encapsulated
   */
  virtual std::shared_ptr<ConduitCapsule> GetChild(std::size_t iChild) = 0;

  /**
   * @brief Get the name of the underlying node
   * @returns A string with the name of the node
   */
  virtual std::string GetName() const = 0;

  ///@{
  /**
   * @brief Type checkers for the data held by the node
   * @returns true if the data type corresponds to the call
   */
  virtual bool IsList() const = 0;
  virtual bool IsInt8() const = 0;
  virtual bool IsInt16() const = 0;
  virtual bool IsInt32() const = 0;
  virtual bool IsInt64() const = 0;
  virtual bool IsUInt8() const = 0;
  virtual bool IsUInt16() const = 0;
  virtual bool IsUInt32() const = 0;
  virtual bool IsUInt64() const = 0;
  virtual bool IsFloat32() const = 0;
  virtual bool IsFloat64() const = 0;
  virtual bool IsString() const = 0;
  ///@}

  ///@{
  /**
   * @brief Check if the data held by the node is external
   */
  virtual bool IsDataExternal() const = 0;
  /**
   * @brief Get a void* pointer to the external data when applicable
   * @returns pointer to external data when relevant and nullptr if not
   */
  virtual void* GetExternalDataPointer() const = 0;
  /**
   * @brief Get the offset with which to read the external data
   * @retuns the offset meta-data
   */
  virtual long GetExternalDataOffset() const = 0;
  /**
   * @brief Get the stride with which to read the external data
   * @retuns the stride meta-data
   */
  virtual long GetExternalDataStride() const = 0;
  /**
   * @brief Get the number of bytes one element of external data occupies
   * @retuns the element size in bytes meta-data
   */
  virtual long GetExternalDataElementBytes() const = 0;
  /**
   * @brief Get the endianness with which to read the external data
   * @retuns the endianness meta-data
   */
  virtual long GetExternalDataEndianness() const = 0;
  ///@}


  ///@{
  /**
   * Get the values of data for conventional data types
   */
  virtual int8_t GetInt8Value() const = 0;
  virtual int16_t GetInt16Value()  const = 0;
  virtual int32_t GetInt32Value()  const = 0;
  virtual int64_t GetInt64Value()  const = 0;
  virtual uint8_t GetUInt8Value()   const = 0;
  virtual uint16_t GetUInt16Value()  const = 0;
  virtual uint32_t GetUInt32Value()  const = 0;
  virtual uint64_t GetUInt64Value()  const = 0;
  virtual float GetFloat32Value() const = 0;
  virtual double GetFloat64Value() const = 0;
  virtual std::string GetStringValue() const = 0;
  ///@}

  ///@{
  /**
   * Get the value arrays of data for conventional data types
   */
  virtual std::vector<int8_t> GetInt8Array() const = 0;
  virtual std::vector<int16_t> GetInt16Array() const = 0;
  virtual std::vector<int32_t> GetInt32Array() const = 0;
  virtual std::vector<int64_t> GetInt64Array() const = 0;
  virtual std::vector<uint8_t> GetUInt8Array()  const = 0;
  virtual std::vector<uint16_t> GetUInt16Array() const = 0;
  virtual std::vector<uint32_t> GetUInt32Array() const = 0;
  virtual std::vector<uint64_t> GetUInt64Array() const = 0;
  virtual std::vector<float> GetFloat32Array() const = 0;
  virtual std::vector<double> GetFloat64Array() const = 0;
  ///@}

};

}

#endif
