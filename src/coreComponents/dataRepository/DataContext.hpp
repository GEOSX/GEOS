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
 * @file DataContext.hpp
 */

#ifndef GEOS_DATAREPOSITORY_DATACONTEXT_HPP_
#define GEOS_DATAREPOSITORY_DATACONTEXT_HPP_

#include "common/DataTypes.hpp"
#include "common/Logger.hpp"
#include "xmlWrapper.hpp"
#include "common/Format.hpp"

namespace geos
{
namespace dataRepository
{


/**
 * @class DataContext
 *
 * DataContext is an abstract class storing contextual information on an object:
 * - Either its line position in a file (if applicable, implementation in DataFileContext),
 * - or its location in the data hierarchy (implementation in GroupContext and WrapperContext).
 * Typically, the target object contains an unique_ptr< DataContext > instance of this class.
 */
class DataContext
{
public:

  /**
   * @brief Construct a new DataContext object.
   * @param targetName the target object name
   */
  DataContext( string const & targetName );

  /**
   * @brief Destroy the DataContext object
   */
  virtual ~DataContext() {}

  /**
   * @return A string that mention all the known informations to retrieve from where the target
   * object comes from.
   */
  virtual string toString() const = 0;

  /**
   * @return Get the target object name
   */
  string getTargetName() const
  { return m_targetName; }
  /**
   * @brief Insert contextual information in the provided stream.
   */
  friend std::ostream & operator<<( std::ostream & os, const DataContext & ctx );

protected:
  // GroupContext & WrapperContext are friend class to be able to access to the protected method on other instances.
  friend class GroupContext;
  friend class WrapperContext;

  /// @see getObjectName()
  string const m_targetName;

  /// This struct exposes the raw data of a DataContext instance that toString() need in order to format it.
  /// This struct lifetime depends on that of the source DataContext. The DataContext is considered constant.
  struct ToStringInfo
  {
    /// the targetName of the DataContext
    string m_targetName;
    /// the file path of the DataFileContext, if it exists (an empty string otherwise)
    string m_filePath;
    /// the file line of the DataFileContext, if it exists (an empty string otherwise)
    size_t m_line = xmlWrapper::xmlDocument::npos;

    /**
     * @brief Construct a new ToStringInfo object from a DataContext that has input file info.
     * @param targetName the target name
     * @param filePath the input file path where the target is declared.
     * @param line the line in the file where the target is declared.
     */
    ToStringInfo( string const & targetName, string const & filePath, size_t line );
    /**
     * @brief Construct a new ToStringInfo object from a DataContext that has no input file info.
     * @param targetName the target name.
     */
    ToStringInfo( string const & targetName );
    /**
     * @return true if a location has been found to declare the target in an input file.
     */
    bool hasInputFileInfo() const
    { return !m_filePath.empty() && m_line != xmlWrapper::xmlDocument::npos; }
  };

  /**
   * @brief This method exposes the raw data of a DataContext, in order to access and format it
   * (notably in toString() implementations that need to access other DataContext instances).
   * @return a ToStringInfo struct that contains the raw data contained in this DataContext instance.
   */
  virtual ToStringInfo getToStringInfo() const = 0;

};

/**
 * @class DataFileContext
 *
 * Stores information to retrieve where a target object has been declared in the input source
 * file (e.g. XML).
 */
class DataFileContext final : public DataContext
{
public:

  /**
   * @brief Construct the file context of a Group from an xml node.
   * @param targetNode the target object xml node
   * @param nodePos the target object xml node position
   */
  DataFileContext( xmlWrapper::xmlNode const & targetNode, xmlWrapper::xmlNodePos const & nodePos );
  /**
   * @brief Construct the file context of a Group from an xml node.
   * @param targetNode the xml node containing the xml attribute
   * @param att the target object xml attribute
   * @param attPos the target object xml attribute position
   */
  DataFileContext( xmlWrapper::xmlNode const & targetNode, xmlWrapper::xmlAttribute const & att,
                   xmlWrapper::xmlAttributePos const & attPos );

  /**
   * @return the target object name followed by the the file and line declaring it.
   */
  string toString() const override;

  /**
   * @return the type name in the source file (XML node tag name / attribute name).
   */
  string getTypeName() const
  { return m_typeName; }

  /**
   * @return the source file path where the target object has been declared.
   */
  string getFilePath() const
  { return m_filePath; }

  /**
   * @return the line (starting from 1) where the target object has been declared in the source file.
   */
  size_t getLine() const
  { return m_line; }

  /**
   * @return the character offset in the line (starting from 1) where the target object has been
   * declared in the source file.
   */
  size_t getOffsetInLine() const
  { return m_offsetInLine; }

  /**
   * @return the character offset (starting from 0) in the source file path where the target object has been
   * declared.
   */
  size_t getOffset() const
  { return m_offset; }

private:

  /// @see getTypeName()
  string const m_typeName;
  /// @see getFilePath()
  string const m_filePath;
  /// @see getLine()
  size_t const m_line;
  /// @see getLineOffset()
  size_t const m_offsetInLine;
  /// @see getOffset()
  size_t const m_offset;

  /**
   * @copydoc DataContext::getToStringInfo()
   */
  ToStringInfo getToStringInfo() const override;

};


} /* namespace dataRepository */
} /* namespace geos */



/**
 * @brief Formatter to be able to directly use a DataContext as a GEOS_FMT() argument.
 * Inherits from formatter<std::string> to reuse its parse() method.
 */
template<>
struct GEOS_FMT_NS::formatter< geos::dataRepository::DataContext > : GEOS_FMT_NS::formatter< std::string >
{
  /**
   * @brief Format the specified DataContext to a string.
   * @param dataContext the DataContext object to format
   * @param ctx formatting state consisting of the formatting arguments and the output iterator
   * @return iterator to the output buffer
   */
  auto format( geos::dataRepository::DataContext const & dataContext, format_context & ctx ) const
  {
    return GEOS_FMT_NS::formatter< std::string >::format( dataContext.toString(), ctx );
  }
};

#endif /* GEOS_DATAREPOSITORY_DATACONTEXT_HPP_ */
