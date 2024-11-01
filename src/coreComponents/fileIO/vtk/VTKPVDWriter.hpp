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

#ifndef GEOS_FILEIO_VTK_VTKPVDWRITER_HPP_
#define GEOS_FILEIO_VTK_VTKPVDWRITER_HPP_

#include "dataRepository/xmlWrapper.hpp"

namespace geos
{
namespace vtk
{

/*!
 * @brief VTK PVD Writer class.
 * @details the PVD file is the root file, it contains the
 * reference to all the VTM files (one VTM file per output time-step).
 * This file is meant to be open with paraview.
 */
class VTKPVDWriter
{
public:
  /*!
   * @brief Constructor
   * @param[in] fileName the file name (with extension)
   */
  explicit VTKPVDWriter( string fileName );

  /*!
   * @brief Set the output file name
   * @param fileName the file name (with extension)
   */
  void setFileName( string fileName );

  /*!
   * @brief Triggers the file output
   */
  void save() const;

  /*!
   * @brief Add a dataset associated to a time-step
   * @param[in] time the time step
   * @param[in] filePath path to the file associated with the time-step
   */
  void addData( real64 time, string const & filePath ) const;

  /*!
   * @brief Removes all the datasets accumulated
   */
  void reinitData() const;

private:

  /// PVD XML file
  xmlWrapper::xmlDocument m_pvdFile;

  /// Name of the XML File
  string m_fileName;
};
}
}

#endif
