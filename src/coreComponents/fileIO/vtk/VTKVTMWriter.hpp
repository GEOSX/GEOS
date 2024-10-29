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

#ifndef GEOS_FILEIO_VTK_VTKVTMWRITER_HPP_
#define GEOS_FILEIO_VTK_VTKVTMWRITER_HPP_

#include "dataRepository/xmlWrapper.hpp"

namespace geos
{
namespace vtk
{

/*!
 * @brief VTM Writer class.
 * @details a VTM file is the root file for one time step. It will contain
 * path to all the component of the mesh (Surfaces, Volumes, Wells etc.).
 */
class VTKVTMWriter
{
public:
  /*!
   * @brief Build the VTM Writer
   * @param[in] filePath path to the file
   */
  explicit VTKVTMWriter( string filePath );

  /*!
   * @brief Triggers the file output
   */
  void write() const;

  /*!
   * @brief Add a dataset block to the VTM file
   * @param[in] blockPath path consisting of intermediate block names that will be created on demand
   * @param[in] dataSetName name of the dataset (leaf level in the multi-block tree)
   * @param[in] filePath path to the dataset file
   */
  void addDataSet( std::vector< string > const & blockPath,
                   string const & dataSetName,
                   string const & filePath ) const;

private:

  /// VTM XML file
  xmlWrapper::xmlDocument m_document;

  /// Handle to the block root node
  xmlWrapper::xmlNode m_blockRoot;

  /// Path to the XML File
  string m_filePath;
};

} // namespace vtk
} // namespace geos

#endif
