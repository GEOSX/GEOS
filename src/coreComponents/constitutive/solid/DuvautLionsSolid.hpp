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
 * @file DuvautLionsSolid.hpp
 * @brief This class implements Duvaut-Lions viscoplasticity model
 *
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_DUVAUTLIONSSOLID_HPP_
#define GEOS_CONSTITUTIVE_SOLID_DUVAUTLIONSSOLID_HPP_

#include "constitutive/solid/SolidBase.hpp"
#include "ElasticIsotropic.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "SolidModelDiscretizationOpsIsotropic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{

// DAMAGE MODEL UPDATES
//
// NOTE: This model uses the m_newStress array to represent the stress in a
//       "non-viscous" elasto-plastic solid.  We then scale the results
//       by the damage factor whenever the true stress is requested through an update
//       function.  The developer should be very cautious if accessing the stress
//       directly through an arrayView, as it does not represent the true stress.
//

template< typename UPDATE_BASE >
class DuvautLionsSolidUpdates : public UPDATE_BASE
{
public:
  template< typename ... PARAMS >
  DuvautLionsSolidUpdates( real64 const & relaxationTime,
                           PARAMS && ... baseParams ):
    UPDATE_BASE( std::forward< PARAMS >( baseParams )... ),
    m_relaxationTime( relaxationTime )
  {}

  using DiscretizationOps = typename UPDATE_BASE::DiscretizationOps; // TODO: typo in anistropic (fix in DiscOps PR)

  using UPDATE_BASE::smallStrainUpdate;
  //using UPDATE_BASE::smallStrainUpdate_ElasticOnly;
  using UPDATE_BASE::saveConvergedState;

  using UPDATE_BASE::viscousStateUpdate;

  using UPDATE_BASE::m_disableInelasticity;

  GEOS_HOST_DEVICE
  //inline
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & timeIncrement,
                          real64 const ( &strainIncrement )[6],
                          real64 ( & stress )[6],
                          real64 ( & stiffness )[6][6] ) const
  {
    real64 trialStress[6];   // Trial stress (elastic predictor)
    real64 elasticStiffness[6][6];  //Elastic stiffness
    real64 timeRatio = 1.0 / (1.0 + timeIncrement / m_relaxationTime);

    for( localIndex i=0; i<6; ++i )
    {
      trialStress[i] = stress[i];
    }

    if( m_disableInelasticity )
    {
      return;
    }


    UPDATE_BASE::smallStrainUpdate_ElasticOnly( k, q, timeIncrement, strainIncrement, trialStress, elasticStiffness );


    UPDATE_BASE::smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness );


    for( localIndex i=0; i<6; ++i )
    {
      stress[i] = timeRatio *  trialStress[i] + (1-timeRatio) * stress[i];
      for( localIndex j=0; j<6; ++j )
      {
        stiffness[i][j] = timeRatio * elasticStiffness[i][j]  + (1 - timeRatio) * stiffness[i][j];
      }
    }

    UPDATE_BASE::saveStress( k, q, stress );
    UPDATE_BASE::viscousStateUpdate( k, q, timeRatio );
    return;
  }

//TODO: modify implementation of smallStrainUpdate to use optimized stiffness -
// this implementation uses full stiffness tensor
//  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  DiscretizationOps & stiffness ) const override final
  {
    this->smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness.m_c );
  }

  real64 const m_relaxationTime;

};



class DuvautLionsBase : public SolidBase
{};

/**
 * @brief Class to represent a rate-dependent Duvaut-Lions material coupled with a plasticity model.
 * It is used as an interface to access all constitutive models relative to the material properties.
 *
 * @tparam SOLID_TYPE type of solid model
 */

template< typename BASE >
class DuvautLionsSolid : public BASE
{
public:

  /// @typedef Alias for DuvautLionsSolidUpdates
  using KernelWrapper = DuvautLionsSolidUpdates< typename BASE::KernelWrapper >;

  DuvautLionsSolid( string const & name, dataRepository::Group * const parent );
  virtual ~DuvautLionsSolid() override = default;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */

  static string catalogName() { return string( "Visco" ) + BASE::m_catalogNameString; }
  virtual string getCatalogName() const override { return catalogName(); }

  real64 relaxationTime() const { return m_relaxationTime; }

  virtual void postInputInitialization() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;


  KernelWrapper createKernelUpdates() const
  {
    return BASE::template createDerivedKernelUpdates< KernelWrapper >( m_relaxationTime );
  }

  struct viewKeyStruct : public BASE::viewKeyStruct
  {
    /// string/key for relaxation time
    static constexpr char const * relaxationTimeString() { return "relaxationTime"; }
  };


protected:
  real64 m_relaxationTime;
};

}
} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_DUVAUTLIONSSOLID_HPP_ */
