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
 * @file CompressibleSolid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_COMPRESSIBLESOLILD_HPP_
#define GEOS_CONSTITUTIVE_SOLID_COMPRESSIBLESOLILD_HPP_

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
class CompressibleSolidUpdates : public CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >
{
public:

  /**
   * @brief Constructor
   */
  CompressibleSolidUpdates( NullModel const & solidModel,
                            PORO_TYPE const & porosityModel,
                            PERM_TYPE const & permModel ):
    CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >( solidModel, porosityModel, permModel )
  {}

  GEOS_HOST_DEVICE
  virtual void updateStateFromPressureAndTemperature( localIndex const k,
                                                      localIndex const q,
                                                      real64 const & pressure,
                                                      real64 const & GEOS_UNUSED_PARAM( pressure_k ),
                                                      real64 const & GEOS_UNUSED_PARAM( pressure_n ),
                                                      real64 const & temperature,
                                                      real64 const & GEOS_UNUSED_PARAM( temperature_k ),
                                                      real64 const & GEOS_UNUSED_PARAM( temperature_n ) ) const override final
  {
    m_porosityUpdate.updateFromPressureAndTemperature( k, q, pressure, temperature );
    real64 const porosity = m_porosityUpdate.getPorosity( k, q );
    m_permUpdate.updateFromPressureAndPorosity( k, q, pressure, porosity );
  }

  GEOS_HOST_DEVICE
  void updateStateFromPressureAndAperture( localIndex const k,
                                           localIndex const q,
                                           real64 const & pressure,
                                           real64 const & oldHydraulicAperture,
                                           real64 const & newHydraulicAperture ) const
  {
    real64 const temperature = 0;
    real64 const dHydraulicAperture_dNormalJump = 1.0;
    m_porosityUpdate.updateFromPressureAndTemperature( k, q, pressure, temperature );
    m_permUpdate.updateFromAperture( k, q, oldHydraulicAperture, newHydraulicAperture, dHydraulicAperture_dNormalJump );
  }

  GEOS_HOST_DEVICE
  void updateStateFromPressureApertureJumpAndTraction( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & pressure,
                                                       real64 const & oldHydraulicAperture,
                                                       real64 const & newHydraulicAperture,
                                                       real64 const & dHydraulicAperture_dNormalJump,
                                                       real64 const ( &dispJump )[3],
                                                       real64 const ( &traction )[3] ) const
  {
    m_porosityUpdate.updateFromPressureAndTemperature( k, q, pressure, 0.0 );
    m_permUpdate.updateFromApertureAndShearDisplacement( k, q, oldHydraulicAperture, newHydraulicAperture, dHydraulicAperture_dNormalJump, pressure, dispJump, traction );
  }

private:
  using CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >::m_solidUpdate;
  using CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >::m_porosityUpdate;
  using CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >::m_permUpdate;

};


/**
 * @brief CompressibleSolidBase class used for dispatch of all Compressible solids.
 */
class CompressibleSolidBase
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
class CompressibleSolid : public CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >
{
public:


  /// Alias for ElasticIsotropicUpdates
  using KernelWrapper = CompressibleSolidUpdates< PORO_TYPE, PERM_TYPE >;

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  CompressibleSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~CompressibleSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "CompressibleSolid" ) + PERM_TYPE::catalogName(); }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }


  /**
   * @brief Create a instantiation of the CompressibleSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of CompressibleSolidUpdates.
   */
  KernelWrapper createKernelUpdates() const
  {

    return CompressibleSolidUpdates< PORO_TYPE, PERM_TYPE >( getSolidModel(),
                                                             getPorosityModel(),
                                                             getPermModel() );
  }
private:
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getSolidModel;
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getPorosityModel;
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getPermModel;

};

}
} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_POROELASTIC_HPP_ */
