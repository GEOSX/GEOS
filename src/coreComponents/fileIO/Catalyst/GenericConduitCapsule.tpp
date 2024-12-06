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
 * @file GenericConduitCapsule.tpp
 */

// source includes
#include "GenericConduitCapsule.hpp"

// external includes
#include <algorithm>

namespace geos
{

template<typename ConduitNodeT>
std::size_t GenericConduitCapsule<ConduitNodeT>::GetNumberOfChildren() const
{
  return this->Node->number_of_children();
}

template<typename ConduitNodeT>
std::shared_ptr<ConduitCapsule> GenericConduitCapsule<ConduitNodeT>::GetChild(std::size_t iChild)
{
  if (iChild >= this->GetNumberOfChildren())
  {
    return nullptr;
  }
  return std::make_shared<GenericConduitCapsule<ConduitNodeT>>(this->Node->child_ptr(iChild));
}

template<typename ConduitNodeT>
std::string GenericConduitCapsule<ConduitNodeT>::GetName() const
{
  return this->Node->name();
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsList() const
{
  return this->Node->dtype().number_of_elements() > 1;
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsInt8() const
{
  return this->Node->dtype().is_int8();
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsInt16() const
{
  return this->Node->dtype().is_int16();
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsInt32() const
{
  return this->Node->dtype().is_int32();
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsInt64() const
{
  return this->Node->dtype().is_int64();
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsUInt8() const
{
  return this->Node->dtype().is_uint8();
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsUInt16() const
{
  return this->Node->dtype().is_uint16();
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsUInt32() const
{
  return this->Node->dtype().is_uint32();
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsUInt64() const
{
  return this->Node->dtype().is_uint64();
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsFloat32() const
{
  return this->Node->dtype().is_float32();
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsFloat64() const
{
  return this->Node->dtype().is_float64();
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsString() const
{
  return this->Node->dtype().is_string();
}

template<typename ConduitNodeT>
bool GenericConduitCapsule<ConduitNodeT>::IsDataExternal() const
{
  return this->Node->is_data_external();
}

template<typename ConduitNodeT>
void* GenericConduitCapsule<ConduitNodeT>::GetExternalDataPointer() const
{
  return this->Node->data_ptr();
}

template<typename ConduitNodeT>
long GenericConduitCapsule<ConduitNodeT>::GetExternalDataOffset() const
{
  return this->Node->dtype().offset();
}

template<typename ConduitNodeT>
long GenericConduitCapsule<ConduitNodeT>::GetExternalDataStride() const
{
  return this->Node->dtype().stride();
}

template<typename ConduitNodeT>
long GenericConduitCapsule<ConduitNodeT>::GetExternalDataElementBytes() const
{
  return this->Node->dtype().element_bytes();
}

template<typename ConduitNodeT>
long GenericConduitCapsule<ConduitNodeT>::GetExternalDataEndianness() const
{
  return this->Node->dtype().endianness();
}

template<typename ConduitNodeT>
int8_t GenericConduitCapsule<ConduitNodeT>::GetInt8Value() const
{
  return this->Node->as_int8();
}

template<typename ConduitNodeT>
int16_t GenericConduitCapsule<ConduitNodeT>::GetInt16Value()  const
{
  return this->Node->as_int16();
}

template<typename ConduitNodeT>
int32_t GenericConduitCapsule<ConduitNodeT>::GetInt32Value()  const
{
  return this->Node->as_int32();
}

template<typename ConduitNodeT>
int64_t GenericConduitCapsule<ConduitNodeT>::GetInt64Value()  const
{
  return this->Node->as_int64();
}

template<typename ConduitNodeT>
uint8_t GenericConduitCapsule<ConduitNodeT>::GetUInt8Value()   const
{
  return this->Node->as_uint8();
}

template<typename ConduitNodeT>
uint16_t GenericConduitCapsule<ConduitNodeT>::GetUInt16Value()  const
{
  return this->Node->as_uint16();
}

template<typename ConduitNodeT>
uint32_t GenericConduitCapsule<ConduitNodeT>::GetUInt32Value()  const
{
  return this->Node->as_uint32();
}

template<typename ConduitNodeT>
uint64_t GenericConduitCapsule<ConduitNodeT>::GetUInt64Value()  const
{
  return this->Node->as_uint64();
}

template<typename ConduitNodeT>
float GenericConduitCapsule<ConduitNodeT>::GetFloat32Value() const
{
  return this->Node->as_float32();
}

template<typename ConduitNodeT>
double GenericConduitCapsule<ConduitNodeT>::GetFloat64Value() const
{
  return this->Node->as_float64();
}

template<typename ConduitNodeT>
std::string GenericConduitCapsule<ConduitNodeT>::GetStringValue() const
{
  return this->Node->as_char8_str();
}

template<typename ConduitNodeT>
std::vector<int8_t> GenericConduitCapsule<ConduitNodeT>::GetInt8Array() const
{
  auto array = this->Node->as_int8_array();
  std::vector<int8_t> v(array.number_of_elements());
  std::generate(v.begin(), v.end(), [i=-1, array] () mutable {i++;return array.element(i);});
  return v;
}

template<typename ConduitNodeT>
std::vector<int16_t> GenericConduitCapsule<ConduitNodeT>::GetInt16Array() const
{
  auto array = this->Node->as_int16_array();
  std::vector<int16_t> v(array.number_of_elements());
  std::generate(v.begin(), v.end(), [i=-1, array] () mutable {i++;return array.element(i);});
  return v;
}

template<typename ConduitNodeT>
std::vector<int32_t> GenericConduitCapsule<ConduitNodeT>::GetInt32Array() const
{
  auto array = this->Node->as_int32_array();
  std::vector<int32_t> v(array.number_of_elements());
  std::generate(v.begin(), v.end(), [i=-1, array] () mutable {i++;return array.element(i);});
  return v;
}

template<typename ConduitNodeT>
std::vector<int64_t> GenericConduitCapsule<ConduitNodeT>::GetInt64Array() const
{
  auto array = this->Node->as_int64_array();
  std::vector<int64_t> v(array.number_of_elements());
  std::generate(v.begin(), v.end(), [i=-1, array] () mutable {i++;return array.element(i);});
  return v;
}

template<typename ConduitNodeT>
std::vector<uint8_t> GenericConduitCapsule<ConduitNodeT>::GetUInt8Array()  const
{
  auto array = this->Node->as_uint8_array();
  std::vector<uint8_t> v(array.number_of_elements());
  std::generate(v.begin(), v.end(), [i=-1, array] () mutable {i++;return array.element(i);});
  return v;
}

template<typename ConduitNodeT>
std::vector<uint16_t> GenericConduitCapsule<ConduitNodeT>::GetUInt16Array() const
{
  auto array = this->Node->as_uint16_array();
  std::vector<uint16_t> v(array.number_of_elements());
  std::generate(v.begin(), v.end(), [i=-1, array] () mutable {i++;return array.element(i);});
  return v;
}

template<typename ConduitNodeT>
std::vector<uint32_t> GenericConduitCapsule<ConduitNodeT>::GetUInt32Array() const
{
  auto array = this->Node->as_uint32_array();
  std::vector<uint32_t> v(array.number_of_elements());
  std::generate(v.begin(), v.end(), [i=-1, array] () mutable {i++;return array.element(i);});
  return v;
}

template<typename ConduitNodeT>
std::vector<uint64_t> GenericConduitCapsule<ConduitNodeT>::GetUInt64Array() const
{
  auto array = this->Node->as_uint64_array();
  std::vector<uint64_t> v(array.number_of_elements());
  std::generate(v.begin(), v.end(), [i=-1, array] () mutable {i++;return array.element(i);});
  return v;
}

template<typename ConduitNodeT>
std::vector<float> GenericConduitCapsule<ConduitNodeT>::GetFloat32Array() const
{
  auto array = this->Node->as_float32_array();
  std::vector<float> v(array.number_of_elements());
  std::generate(v.begin(), v.end(), [i=-1, array] () mutable {i++;return array.element(i);});
  return v;
}

template<typename ConduitNodeT>
std::vector<double> GenericConduitCapsule<ConduitNodeT>::GetFloat64Array() const
{
  auto array = this->Node->as_float64_array();
  std::vector<double> v(array.number_of_elements());
  std::generate(v.begin(), v.end(), [i=-1, array] () mutable {i++;return array.element(i);});
  return v;
}

}
