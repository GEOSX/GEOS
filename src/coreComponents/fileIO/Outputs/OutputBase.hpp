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
 * @file OutputBase.hpp
 */
#ifndef GEOS_FILEIO_OUTPUTS_OUTPUTBASE_HPP_
#define GEOS_FILEIO_OUTPUTS_OUTPUTBASE_HPP_

#include "dataRepository/Group.hpp"
#include "dataRepository/ExecutableGroup.hpp"


namespace geos
{

/**
 * @class OutputBase
 *
 * A base class for output types.
 */
class OutputBase : public ExecutableGroup
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  explicit OutputBase( string const & name, Group * const parent );

  /// Destructor
  virtual ~OutputBase() override;

  /**
   * @brief Setter for the output directory
   * @param  outputDir The output directory
   **/
  static void setOutputDirectory( string const & outputDir );

  /**
   * @brief Getter for the output directory
   * @return The output directory
   **/
  static string getOutputDirectory() {return m_outputDirectory;}

  /**
   * @brief Setter for the file name root
   * @param root The file name root
   **/
  static void setFileNameRoot( string const & root );

  /**
   * @brief Getter for the file name root
   * @return The file name root
   **/
  static string getFileNameRoot() { return m_fileNameRoot; }

  /// Method for setting up output directories.
  virtual void setupDirectoryStructure();

  // Catalog interface
  /// @cond DO_NOT_DOCUMENT
  using CatalogInterface = dataRepository::CatalogInterface< OutputBase, string const &, Group * const >;
  static CatalogInterface::CatalogType & getCatalog();

  // Catalog view keys
  struct viewKeysStruct
  {
    static constexpr auto childDirectoryString = "childDirectory";
    static constexpr auto parallelThreadsString = "parallelThreads";
  } outputBaseViewKeys;
  /// @endcond

  /**
   * @brief Get the path of the child directory where output will be written
   * @return The directory path
   **/
  string childDirectory() const { return m_childDirectory; }

  /**
   * @brief Get the number of parallel threads to use to write plotfiles
   * @return The number of threads
   **/
  integer parallelThreads() const { return m_parallelThreads; }

protected:
  /**
   * @brief Do initialization prior to calling initialization operations
   *        on the subgroups.
   * @param group The root group
   **/
  virtual void initializePreSubGroups() override;

private:
  string m_childDirectory;
  integer m_parallelThreads;

  static string m_outputDirectory;
  static string m_fileNameRoot;

};


} /* namespace geos */

#endif /* GEOS_FILEIO_OUTPUTS_OUTPUTBASE_HPP_ */
