/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ExternalMeshGeneratorBase.hpp
 */

#ifndef GEOS_MESH_GENERATORS_EXTERNALMESHGENERATORBASE_HPP
#define GEOS_MESH_GENERATORS_EXTERNALMESHGENERATORBASE_HPP

#include "mesh/generators/MeshGeneratorBase.hpp"

namespace geos
{

/**
 * @brief Base class for external mesh generators (importers).
 */
class ExternalMeshGeneratorBase : public MeshGeneratorBase
{
public:

  /**
   * @brief Constructor.
   * @param[in] name name of the object
   * @param[in] parent the parent Group pointer
   */
  ExternalMeshGeneratorBase( const string & name,
                             Group * const parent );

protected:

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * filePathString() { return "file"; }
    constexpr static char const * scaleString() { return "scale"; }
    constexpr static char const * translateString() { return "translate"; }
    constexpr static char const * volumicFieldsToImportString() { return "fieldsToImport"; }
    constexpr static char const * volumicFieldsInGEOSXString() { return "fieldNamesInGEOSX"; }
    constexpr static char const * surfacicFieldsToImportString() { return "surfacicFieldsToImport"; }
    constexpr static char const * surfacicFieldsInGEOSXString() { return "surfacicFieldsInGEOSX"; }
  };
  /// @endcond

  void postInputInitialization() override;

  /// Path to the mesh file
  Path m_filePath;

  /// Translation vector that will be applied to the point coordinates (prior to scaling)
  R1Tensor m_translate;

  /// Scale factor that will be applied to the point coordinates (after translation)
  R1Tensor m_scale;

  /// Names of the fields to be copied from an external reader into GEOSX data structure
  array1d< string > m_volumicFieldsToImport;

  /// String array of the GEOSX user declared volumic fields
  array1d< string > m_volumicFieldsInGEOSX;

  /// Names of the surfacic fields to be copied from an external reader into GEOSX data structure
  array1d< string > m_surfacicFieldsToImport;

  /// String array of the GEOSX user declared surfacic fields
  array1d< string > m_surfacicFieldsInGEOSX;
};

} // namespace geos

#endif //GEOS_MESH_GENERATORS_EXTERNALMESHGENERATORBASE_HPP
