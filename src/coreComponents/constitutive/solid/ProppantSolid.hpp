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
 * @file ProppantSolid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_PROPPANTSOLID_HPP_
#define GEOS_CONSTITUTIVE_SOLID_PROPPANTSOLID_HPP_

#include "constitutive/solid/CoupledSolid.hpp"
#include "constitutive/NullModel.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 *
 *
 * @tparam PORO_TYPE type of the porosity model
 * @tparam PERM_TYPE type of the permeability model
 */
template< typename PORO_TYPE,
          typename PERM_TYPE >
class ProppantSolidUpdates : public CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >
{
public:

  /**
   * @brief Constructor
   */
  ProppantSolidUpdates( NullModel const & solidModel,
                        PORO_TYPE const & porosityModel,
                        PERM_TYPE const & permModel ):
    CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >( solidModel, porosityModel, permModel )
  {}

  GEOS_HOST_DEVICE
  void updateStateFromApertureAndProppantVolumeFraction( localIndex const k,
                                                         localIndex const q,
                                                         real64 const & oldHydraulicAperture,
                                                         real64 const & newHydraulicAperture,
                                                         real64 const & proppantPackVolumeFraction ) const
  {
    m_porosityUpdate.updateFromProppantVolumeFraction( k, q, proppantPackVolumeFraction );
    real64 const dHydraulicAperture_dNormalJump = 1.0;
    m_permUpdate.updateFromApertureAndProppantVolumeFraction( k, q, oldHydraulicAperture, newHydraulicAperture, dHydraulicAperture_dNormalJump, proppantPackVolumeFraction );
  }

private:
  using CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >::m_solidUpdate;
  using CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >::m_porosityUpdate;
  using CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >::m_permUpdate;

};


/**
 * @brief ProppantSolidBase class used for dispatch of all Proppant solids.
 */
class ProppantSolidBase
{};


/**
 * @brief Class to represent a porous material for flow simulations.
 * It is used as an interface to access all constitutive models relative to the properties of a porous material
 * for flow only simulations.
 *
 * @tparam PORO_TYPE type of porosity model
 * @tparam PERM_TYPE type of the permeability model
 */

template< typename PORO_TYPE,
          typename PERM_TYPE >
class ProppantSolid : public CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >
{
public:


  /// Alias for ElasticIsotropicUpdates
  using KernelWrapper = ProppantSolidUpdates< PORO_TYPE, PERM_TYPE >;

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  ProppantSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~ProppantSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "ProppantSolid" ) + PERM_TYPE::catalogName(); }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }


  /**
   * @brief Create an instantiation of the ProppantSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of ProppantSolidUpdates.
   */
  KernelWrapper createKernelUpdates() const
  {

    return ProppantSolidUpdates< PORO_TYPE, PERM_TYPE >( getSolidModel(),
                                                         getPorosityModel(),
                                                         getPermModel() );
  }
private:
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getSolidModel;
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getPorosityModel;
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getPermModel;

};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_PROPPANTSOLID_HPP_ */
