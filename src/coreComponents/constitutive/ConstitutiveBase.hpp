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

/**
 * @file ConstitutiveBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_
#define GEOS_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_

#include "dataRepository/ObjectCatalog.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"

namespace geos
{

namespace constitutive
{


class ConstitutiveBase : public dataRepository::Group
{
public:

  /**
   * Single point of reference for generated constitutive field naming convention
   *
   * @param prefix name prefix (e.g. constitutive model name)
   * @param name actual field name
   * @return prefixed field name that is used to access data
   */
  inline static string makeFieldName( string const & prefix, string const & name ) { return prefix + "_" + name; }


  ConstitutiveBase( string const & name,
                    Group * const parent );

  virtual ~ConstitutiveBase() override = default;

  /**
   * @brief create a clone of this constitutive model
   * @param[in]  name   The name of the clone in the repository
   * @param[in]  parent A pointer to the group that contains the instance of the new clone
   * @param[out] clone  A reference to a unique_ptr  that will hold the clone.
   */
  virtual std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                            Group * const parent ) const;


  /// Save state data in preparation for next timestep
  virtual void saveConvergedState() const
  {}

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// @typedef An alias for the ConstitutiveBase catalog
  using CatalogInterface = dataRepository::CatalogInterface< ConstitutiveBase, string const &, Group * const >;

  /**
   * @brief Singleton accessor for catalog
   * @return
   */
  static typename CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief function to return the catalog name of the derived class
   * @return a string that contains the catalog name of the derived class
   */
  virtual string getCatalogName() const = 0;

  /**
   * @brief Get full name of the model.
   * @return full name, consisting of XML (catalog) name and actual model name, and if possible
   * followed by the xml file and line.
   */
  string getFullName() const
  {
    return getCatalogName() + " " + getDataContext().toString();
  }

  ///@}

  /**
   * @brief Allocate constitutive data and make views to data on parent objects
   * @param[in] parent reference to the group that holds the constitutive relation
   * @param[in] numConstitutivePointsPerParentIndex number of quadrature points
   *
   * This function does 2 things:
   *   1) Allocate data according to the size of parent and numConstitutivePointsPerParentIndex
   *   2) Create wrappers to the constitutive data in the parent for easier access
   */
  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex );

  struct viewKeyStruct
  {};

  localIndex numQuadraturePoints() const { return m_numQuadraturePoints; }

  virtual std::vector< string > getSubRelationNames() const { return {}; }

  /**
   * @brief Helper function to register field on a constitutive model
   * @tparam FIELD_TRAIT the type of field
   * @param[in] fieldTrait the field struct corresponding to the object being registered
   * @param[in] newObject a pointer to the object that is being registered
   * @return A reference to the newly registered/created Wrapper
   * TODO: Remove duplicated code with ObjectManagerBase
   */
  template< typename FIELD_TRAIT >
  dataRepository::Wrapper< typename FIELD_TRAIT::type > & registerField( FIELD_TRAIT const & fieldTrait,
                                                                         typename FIELD_TRAIT::type * newObject )
  {
    return registerWrapper( fieldTrait.key(), newObject ).
             setApplyDefaultValue( fieldTrait.defaultValue() ).
             setPlotLevel( FIELD_TRAIT::plotLevel ).
             setRestartFlags( FIELD_TRAIT::restartFlag ).
             setDescription( FIELD_TRAIT::description );
  }

  /**
   * @brief Get a wrapper associated with a trait from the constitutive model
   * @tparam FIELD_TRAIT The trait that holds the type and key of the field
   *   to be retrieved from this constitutive model
   * @return A const reference to a view to const wrapper.
   * TODO: Remove duplicated code with ObjectManagerBase
   */
  template< typename FIELD_TRAIT >
  dataRepository::Wrapper< typename FIELD_TRAIT::type > const & getField() const
  {
    return this->getWrapper< typename FIELD_TRAIT::type >( FIELD_TRAIT::key() );
  }

  /**
   * @brief Get a wrapper associated with a trait from the constitutive model
   * @tparam FIELD_TRAIT The trait that holds the type and key of the field
   *   to be retrieved from this constitutive model
   * @return A reference to the wrapper.
   * TODO: Remove duplicated code with ObjectManagerBase
   */
  template< typename FIELD_TRAIT >
  dataRepository::Wrapper< typename FIELD_TRAIT::type > & getField()
  {
    return this->getWrapper< typename FIELD_TRAIT::type >( FIELD_TRAIT::key() );
  }

private:

  localIndex m_numQuadraturePoints;
};

}
}
#endif
