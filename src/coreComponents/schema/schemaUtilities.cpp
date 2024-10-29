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
 * @file schemaUtilities.cpp
 */

/// Source includes
#include "schemaUtilities.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/InputFlags.hpp"
#include "dataRepository/Wrapper.hpp"

/// System includes
#include <regex>

namespace geos
{

using namespace dataRepository;

namespace schemaUtilities
{


void ConvertDocumentationToSchema( string const & fname,
                                   Group * const group,
                                   integer documentationType )
{
  GEOS_LOG_RANK_0( "Generating XML Schema..." );

  string schemaBase=
    "<?xml version=\"1.1\" encoding=\"ISO-8859-1\" ?>\
  <xsd:schema xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\">\
  <xsd:annotation>\
  <xsd:documentation xml:lang=\"en\">GEOSX Input Schema</xsd:documentation>\
  </xsd:annotation>\
  </xsd:schema>";

  xmlWrapper::xmlDocument schemaTree;
  schemaTree.loadString( schemaBase );
  xmlWrapper::xmlNode schemaRoot = schemaTree.getChild( "xsd:schema" );

  // Build the simple schema types
  GEOS_LOG_RANK_0( "  Basic datatypes" );
  BuildSimpleSchemaTypes( schemaRoot );

  // Recursively build the schema from the data structure skeleton
  GEOS_LOG_RANK_0( "  Data structure layout" );
  SchemaConstruction( *group, schemaRoot, schemaRoot, documentationType );

  // Write the schema to file
  GEOS_LOG_RANK_0( "  Saving file" );
  schemaTree.saveFile( fname );

  GEOS_LOG_RANK_0( "  Done!" );
}

string getSchemaTypeName( string_view rtTypeName )
{
  string const sanitizedName = std::regex_replace( string( rtTypeName ), std::regex( "::" ), "_" );

  // Note: Some type names involving strings can vary on compiler and be ugly.  Convert these to "string"
  auto constexpr typeMergingRegex = "std_(__cxx11_basic_)?string(<\\s*char,\\s*std_char_traits<char>,\\s*std_allocator<char>\\s*>)?";
  string xmlSafeName = std::regex_replace( sanitizedName,
                                           std::regex( typeMergingRegex ),
                                           "string" );

  // Replace '<', '>', spaces and ',' because the schema does not like them
  xmlSafeName = std::regex_replace( xmlSafeName, std::regex( "<" ), "_lt_" );
  xmlSafeName = std::regex_replace( xmlSafeName, std::regex( ">" ), "_gt_" );
  xmlSafeName = std::regex_replace( xmlSafeName, std::regex( "," ), "_cm_" );
  xmlSafeName = std::regex_replace( xmlSafeName, std::regex( " " ), "-" );

  return xmlSafeName;
}

void AppendSimpleType( xmlWrapper::xmlNode & schemaRoot,
                       string const & name,
                       string const & regex )
{
  string const advanced_match_string = ".*[\\[\\]`$].*|";

  xmlWrapper::xmlNode newNode = schemaRoot.append_child( "xsd:simpleType" );
  newNode.append_attribute( "name" ) = name.c_str();
  xmlWrapper::xmlNode restrictionNode = newNode.append_child( "xsd:restriction" );
  restrictionNode.append_attribute( "base" ) = "xsd:string";
  xmlWrapper::xmlNode patternNode = restrictionNode.append_child( "xsd:pattern" );

  // Handle the default regex
  if( regex.empty() )
  {
    GEOS_WARNING( "schema regex not defined for " << name );
    patternNode.append_attribute( "value" ) = "(?s).*";
  }
  else
  {
    string const patternString = advanced_match_string + regex;
    patternNode.append_attribute( "value" ) = patternString.c_str();
  }
}

void BuildSimpleSchemaTypes( xmlWrapper::xmlNode schemaRoot )
{
  auto const regexes = rtTypes::createBasicTypesRegexMap();
  for( auto const & [typeName, regex] : regexes )
  {
    AppendSimpleType( schemaRoot, getSchemaTypeName( typeName ), regex.m_regexStr );
  }
}


void SchemaConstruction( Group & group,
                         xmlWrapper::xmlNode schemaRoot,
                         xmlWrapper::xmlNode schemaParent,
                         integer documentationType )
{
  // Get schema details
  InputFlags schemaType = group.getInputFlags();

  if((schemaType != InputFlags::INVALID) || (documentationType == 1))
  {
    string const & targetName = group.getName();
    string typeName = targetName + "Type";

    // Insert the schema node if not present, then iterate over children
    if( schemaParent.find_child_by_attribute( "xsd:element", "name", targetName.c_str()).empty())
    {
      // Add the entries to the current and root nodes
      xmlWrapper::xmlNode targetIncludeNode = schemaParent.append_child( "xsd:element" );
      targetIncludeNode.append_attribute( "name" ) = targetName.c_str();
      targetIncludeNode.append_attribute( "type" ) = typeName.c_str();

      // Add occurence conditions
      if((schemaType == InputFlags::REQUIRED_NONUNIQUE) || (schemaType == InputFlags::REQUIRED))
      {
        targetIncludeNode.append_attribute( "minOccurs" ) = "1";
      }
      if((schemaType == InputFlags::OPTIONAL) || (schemaType == InputFlags::REQUIRED))
      {
        targetIncludeNode.append_attribute( "maxOccurs" ) = "1";
      }

      // Insert a new type into the root node if not present
      xmlWrapper::xmlNode targetTypeDefNode = schemaRoot.find_child_by_attribute( "xsd:complexType", "name", typeName.c_str());
      if( targetTypeDefNode.empty())
      {
        targetTypeDefNode = schemaRoot.append_child( "xsd:complexType" );
        targetTypeDefNode.append_attribute( "name" ) = typeName.c_str();
      }

      // Add subgroups
      if( group.numSubGroups() > 0 )
      {
        // Children are defined in a choice node
        xmlWrapper::xmlNode targetChoiceNode = targetTypeDefNode.child( "xsd:choice" );
        if( targetChoiceNode.empty() )
        {
          targetChoiceNode = targetTypeDefNode.prepend_child( "xsd:choice" );
          targetChoiceNode.append_attribute( "minOccurs" ) = "0";
          targetChoiceNode.append_attribute( "maxOccurs" ) = "unbounded";
        }

        // Get a list of the subgroup names in alphabetic order
        // Note: this is necessary because the order that objects
        //       are registered to catalogs may vary by compiler
        std::set< string > subGroupNames;
        for( auto & subGroupPair : group.getSubGroups())
        {
          subGroupNames.insert( subGroupPair.first );
        }

        // Add children of the group
        for( string subName : subGroupNames )
        {
          Group & subGroup = group.getGroup( subName );
          InputFlags subSchemaType = subGroup.getInputFlags();

          if( ( documentationType == 0 ) & (( subSchemaType == InputFlags::REQUIRED_NONUNIQUE ) || ( subSchemaType == InputFlags::OPTIONAL_NONUNIQUE )))
          {
            // Enforce uniqueness of element names
            // Note: this must be done at the parent element level
            xmlWrapper::xmlNode uniqueNameNode = targetIncludeNode.append_child( "xsd:unique" );
            string path = std::regex_replace( group.getPath(), std::regex( "/" ), "" );
            string uniqueNameNodeStr =  path + subName + "UniqueName";
            uniqueNameNode.append_attribute( "name" ) = uniqueNameNodeStr.c_str();
            xmlWrapper::xmlNode uniqueNameSelector = uniqueNameNode.append_child( "xsd:selector" );
            uniqueNameSelector.append_attribute( "xpath" ) = subName.c_str();
            xmlWrapper::xmlNode uniqueNameField = uniqueNameNode.append_child( "xsd:field" );
            uniqueNameField.append_attribute( "xpath" ) = "@name";
          }

          SchemaConstruction( subGroup, schemaRoot, targetChoiceNode, documentationType );
        }
      }

      // Add schema deviations
      group.setSchemaDeviations( schemaRoot, targetTypeDefNode, documentationType );

      // Add attributes
      // Note: wrappers that were added to this group by another group
      //       may end up in different order.  To avoid this, add them
      //       into the schema in alphabetic order.
      std::set< string > groupWrapperNames;
      for( auto & wrapperPair : group.wrappers())
      {
        groupWrapperNames.insert( wrapperPair.first );
      }

      for( string attributeName : groupWrapperNames )
      {
        WrapperBase & wrapper = group.getWrapperBase( attributeName );
        InputFlags flag = wrapper.getInputFlag();

        if(( flag > InputFlags::FALSE ) != ( documentationType == 1 ))
        {
          // Ignore duplicate copies of attributes
          if( targetTypeDefNode.find_child_by_attribute( "xsd:attribute", "name", attributeName.c_str()).empty())
          {
            // Write any additional documentation that isn't expected by the .xsd format in a comment
            // Attribute description
            string const description = wrapper.getDescription();
            string commentString = attributeName + " => ";

            if( !description.empty())
            {
              commentString += description;
            }
            else
            {
              commentString += "(no description available)";
            }

            // List of objects that registered this field
            std::set< string > const & registrars = wrapper.getRegisteringObjects();
            if( !registrars.empty() )
            {
              commentString += " => " + stringutilities::join( registrars.begin(), registrars.end(), ", " );
            }

            xmlWrapper::xmlNode commentNode = targetTypeDefNode.append_child( xmlWrapper::xmlNodeType::node_comment );
            commentNode.set_value( commentString.c_str());


            // Write the valid schema attributes
            // Basic attributes
            xmlWrapper::xmlNode attributeNode = targetTypeDefNode.append_child( "xsd:attribute" );
            attributeNode.append_attribute( "name" ) = attributeName.c_str();

            string const schemaTypeName = getSchemaTypeName( wrapper.getRTTypeName() );
            GEOS_LOG_VAR( schemaTypeName );
            attributeNode.append_attribute( "type" ) = schemaTypeName.c_str();

            // Check if the attribute has a previously unseen non-simple type with a custom validation regex
            if( schemaRoot.find_child_by_attribute( "xsd:simpleType", "name", schemaTypeName.c_str() ).empty() )
            {
              string const & regex = wrapper.getTypeRegex().m_regexStr;
              if( !regex.empty() )
              {
                // Append a new simpleType with a custom regex
                AppendSimpleType( schemaRoot, schemaTypeName, regex );
              }
            }

            // (Optional) Default Value
            if( (flag == InputFlags::OPTIONAL_NONUNIQUE) || (flag == InputFlags::REQUIRED_NONUNIQUE))
            {
              GEOS_LOG_RANK_0( attributeName << " has an invalid input flag" );
              GEOS_ERROR( "SchemaConstruction: duplicate xml attributes are not allowed" );
            }
            else if( flag == InputFlags::OPTIONAL )
            {
              if( wrapper.hasDefaultValue() )
              {
                attributeNode.append_attribute( "default" ) = wrapper.getDefaultValueString().c_str();
              }
            }
            else if( documentationType == 0 )
            {
              attributeNode.append_attribute( "use" ) = "required";
            }
          }
        }
      }

      // Elements that are nonunique require the use of the name attribute
      if(((schemaType == InputFlags::REQUIRED_NONUNIQUE) || (schemaType == InputFlags::OPTIONAL_NONUNIQUE)) && (documentationType == 0))
      {
        // Only add this attribute if not present
        if( targetTypeDefNode.find_child_by_attribute( "xsd:attribute", "name", "name" ).empty())
        {
          xmlWrapper::xmlNode commentNode = targetTypeDefNode.append_child( xmlWrapper::xmlNodeType::node_comment );
          commentNode.set_value( "name => A name is required for any non-unique nodes" );

          xmlWrapper::xmlNode attributeNode = targetTypeDefNode.append_child( "xsd:attribute" );
          string const schemaTypeName = getSchemaTypeName( rtTypes::CustomTypes::groupName );
          attributeNode.append_attribute( "name" ) = "name";
          attributeNode.append_attribute( "type" ) = schemaTypeName.c_str();
          attributeNode.append_attribute( "use" ) = "required";
        }
      }
    }
  }
}

} /// namespace schemaUtilities
} /// namespace geos
