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
 * @file ConstitutivePassThru.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_
#define GEOS_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_

#include "ConstitutivePassThruHandler.hpp"
#include "NullModel.hpp"
#include "ContinuumBase.hpp"
#include "gas/Gas.hpp"
#include "solid/DamageVolDev.hpp"
#include "solid/DamageSpectral.hpp"
#include "solid/DruckerPrager.hpp"
#include "solid/DruckerPragerExtended.hpp"
#include "solid/PerfectlyPlastic.hpp"
#include "solid/ModifiedCamClay.hpp"
#include "solid/DelftEgg.hpp"
#include "solid/DuvautLionsSolid.hpp"
#include "solid/ElasticIsotropic.hpp"
#include "solid/ElasticIsotropicPressureDependent.hpp"
#include "solid/ElasticTransverseIsotropic.hpp"
#include "solid/ElasticTransverseIsotropicPressureDependent.hpp"
#include "solid/Geomechanics.hpp"
#include "solid/Graphite.hpp"
#include "solid/ElasticOrthotropic.hpp"
#include "solid/Hyperelastic.hpp"
#include "solid/HyperelasticMMS.hpp"
#include "solid/PorousSolid.hpp"
#include "solid/CompressibleSolid.hpp"
#include "solid/ProppantSolid.hpp"
#include "solid/StrainHardeningPolymer.hpp"
#include "solid/Chiumenti.hpp"
#include "solid/CeramicDamage.hpp"
#include "solid/VonMisesJ.hpp"
#include "solid/porosity/PressurePorosity.hpp"
#include "solid/porosity/ProppantPorosity.hpp"
#include "permeability/ConstantPermeability.hpp"
#include "permeability/CarmanKozenyPermeability.hpp"
#include "permeability/ExponentialDecayPermeability.hpp"
#include "permeability/ParallelPlatesPermeability.hpp"
#include "permeability/PressurePermeability.hpp"
#include "permeability/ProppantPermeability.hpp"
#include "permeability/SlipDependentPermeability.hpp"
#include "permeability/WillisRichardsPermeability.hpp"


namespace geos
{
namespace constitutive
{

/**
 * @struct ConstitutivePassThru
 * @brief Struct to facilitate launching of lambda functions with a compile
 *   time knowledge of what constitutive model is used.
 *
 * This struct works by implementing an if-else or switch-case block for a
 * specific constitutive base type, and executing the lambda passing it a
 * casted pointer to the constitutive relation.
 */
template< typename BASETYPE >
struct ConstitutivePassThru;

/**
 * Specialization for models that derive from ElasticIsotropic.
 */
template<>
struct ConstitutivePassThru< ElasticIsotropic >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< ElasticIsotropic >::execute( constitutiveRelation,
                                                              std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for models that derive from Hyperelastic.
 */
template<>
struct ConstitutivePassThru< Hyperelastic >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< HyperelasticMMS >::execute( constitutiveRelation,
                                                              std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for models that derive from HyperelasticMMS.
 */
template<>
struct ConstitutivePassThru< HyperelasticMMS >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< HyperelasticMMS >::execute( constitutiveRelation,
                                                              std::forward< LAMBDA >( lambda ) );
  }
};


/**
 * Specialization for models that derive from SolidBase.
 */
template<>
struct ConstitutivePassThru< SolidBase >
{

  // NOTE: The switch order here can be fragile if a model derives from another
  //       model, as the dynamic_cast will also cast to a base version.
  //       Models should be ordered such that children come before parents.
  //       For example, DruckerPrager before ElasticIsotropic, DamageVolDev before
  //       Damage, etc.

  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< DamageSpectral< ElasticIsotropic >,
                                 DamageVolDev< ElasticIsotropic >,
                                 Damage< ElasticIsotropic >,
                                 DuvautLionsSolid< DruckerPrager >,
                                 DuvautLionsSolid< DruckerPragerExtended >,
                                 DuvautLionsSolid< ModifiedCamClay >,
                                 DruckerPragerExtended,
                                 ModifiedCamClay,
                                 DelftEgg,
                                 DruckerPrager,
                                 ElasticIsotropic,
                                 ElasticTransverseIsotropic,
                                 ElasticIsotropicPressureDependent,
                                 ElasticOrthotropic >::execute( constitutiveRelation,
                                                                std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * @struct ConstitutivePassThruMPM
 */
template< typename BASETYPE >
struct ConstitutivePassThruMPM;

/**
 * Specialization for models that derive from ContinuumBase that are used by the MPM solver.
 * NOTE: this is only a temporary dispatch to reduce the compilation time.
 */
template<>
struct ConstitutivePassThruMPM< ContinuumBase >
{

  // NOTE: The switch order here can be fragile if a model derives from another
  //       model, as the dynamic_cast will also cast to a base version.
  //       Models should be ordered such that children come before parents.
  //       For example, DruckerPrager before ElasticIsotropic, DamageVolDev before
  //       Damage, etc.

  template< typename LAMBDA >
  static
  void execute( ContinuumBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< Graphite,
                                 Geomechanics,
                                 CeramicDamage,
                                 Chiumenti,
                                 StrainHardeningPolymer,
                                 PerfectlyPlastic,
                                 ElasticTransverseIsotropicPressureDependent,
                                 ElasticTransverseIsotropic,
                                 VonMisesJ,
                                 ElasticIsotropic,
                                 Hyperelastic,
                                 HyperelasticMMS,
                                 Gas >::execute( constitutiveRelation,
                                                 std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * @struct ConstitutivePassThruTriaxialDriver
 */
template< typename BASETYPE >
struct ConstitutivePassThruTriaxialDriver;

/**
 * Specialization for models that derive from SolidBase.
 * NOTE: this is only a temporary dispatch to reduce the compilation time.
 */
template<>
struct ConstitutivePassThruTriaxialDriver< SolidBase >
{

  // NOTE: The switch order here can be fragile if a model derives from another
  //       model, as the dynamic_cast will also cast to a base version.
  //       Models should be ordered such that children come before parents.
  //       For example, DruckerPrager before ElasticIsotropic, DamageVolDev before
  //       Damage, etc.
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< DamageSpectral< ElasticIsotropic >,
                                 DamageVolDev< ElasticIsotropic >,
                                 Damage< ElasticIsotropic >,
                                 DuvautLionsSolid< DruckerPrager >,
                                 DuvautLionsSolid< DruckerPragerExtended >,
                                 DuvautLionsSolid< ModifiedCamClay >,
                                 DruckerPragerExtended,
                                 ModifiedCamClay,
                                 DelftEgg,
                                 DruckerPrager,
                                 ElasticIsotropic,
                                 ElasticTransverseIsotropic,
                                 ElasticIsotropicPressureDependent,
                                 ElasticOrthotropic >::execute( constitutiveRelation,
                                                                std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for the NullModel.
 */
template<>
struct ConstitutivePassThru< NullModel >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    if( auto * const ptr = dynamic_cast< NullModel * >( &constitutiveRelation ) )
    {
      lambda( *ptr );
    }
    else
    {
      GEOS_ERROR( "ConstitutivePassThru< NullModel >::execute failed on constitutive relation "
                  << constitutiveRelation.getDataContext() << " with type "
                  << LvArray::system::demangleType( constitutiveRelation ) );
    }
  }
};


/**
 * Specialization for the PorousSolid< ElasticIsotropic > model.
 */
template<>
struct ConstitutivePassThru< PorousSolid< ElasticIsotropic > >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    if( auto * const ptr = dynamic_cast< PorousSolid< ElasticIsotropic > * >( &constitutiveRelation ) )
    {
      lambda( *ptr );
    }
    else
    {
      GEOS_ERROR( "ConstitutivePassThru< PorousSolid< ElasticIsotropic > >::execute failed on constitutive relation "
                  << constitutiveRelation.getDataContext() << " with type "
                  << LvArray::system::demangleType( constitutiveRelation ) );
    }
  }
};

/**
 * Specialization for the Damage models.
 */
template<>
struct ConstitutivePassThru< DamageBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation,
                       LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< DamageSpectral< ElasticIsotropic >,
                                 DamageVolDev< ElasticIsotropic >,
                                 Damage< ElasticIsotropic > >::execute( constitutiveRelation,
                                                                        std::forward< LAMBDA >( lambda ) );
  }
};



/**
 * Specialization for the PorousSolid models.
 */
template<>
struct ConstitutivePassThru< PorousSolidBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< PorousSolid< DruckerPragerExtended >,
                                 PorousSolid< ModifiedCamClay >,
                                 PorousSolid< DelftEgg >,
                                 PorousSolid< DruckerPrager >,
                                 PorousSolid< DuvautLionsSolid< DruckerPrager > >,
                                 PorousSolid< DuvautLionsSolid< DruckerPragerExtended > >,
                                 PorousSolid< DuvautLionsSolid< ModifiedCamClay > >,
                                 PorousSolid< ElasticIsotropic >,
                                 PorousSolid< ElasticTransverseIsotropic >,
                                 PorousSolid< ElasticIsotropicPressureDependent >,
                                 PorousSolid< ElasticOrthotropic >,
                                 PorousSolid< DamageSpectral< ElasticIsotropic > >,
                                 PorousSolid< DamageVolDev< ElasticIsotropic > >,
                                 PorousSolid< Damage< ElasticIsotropic > > >::execute( constitutiveRelation,
                                                                                       std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for the CompressibleSolid models.
 */
template<>
struct ConstitutivePassThru< CompressibleSolidBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< CompressibleSolid< PressurePorosity, ConstantPermeability >,
                                 CompressibleSolid< PressurePorosity, CarmanKozenyPermeability >,
                                 CompressibleSolid< PressurePorosity, ExponentialDecayPermeability >,
                                 CompressibleSolid< PressurePorosity, ParallelPlatesPermeability >,
                                 CompressibleSolid< PressurePorosity, PressurePermeability >,
                                 CompressibleSolid< PressurePorosity, SlipDependentPermeability >,
                                 CompressibleSolid< PressurePorosity, WillisRichardsPermeability >
                                 >::execute( constitutiveRelation,
                                             std::forward< LAMBDA >( lambda ) );
  }

  template< typename LAMBDA >
  static void execute( ConstitutiveBase const & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< CompressibleSolid< PressurePorosity, ConstantPermeability >,
                                 CompressibleSolid< PressurePorosity, CarmanKozenyPermeability >,
                                 CompressibleSolid< PressurePorosity, ExponentialDecayPermeability >,
                                 CompressibleSolid< PressurePorosity, ParallelPlatesPermeability >,
                                 CompressibleSolid< PressurePorosity, PressurePermeability >,
                                 CompressibleSolid< PressurePorosity, SlipDependentPermeability >,
                                 CompressibleSolid< PressurePorosity, WillisRichardsPermeability >
                                 >::execute( constitutiveRelation,
                                             std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for the ProppantModel.
 */
template<>
struct ConstitutivePassThru< ProppantSolid< ProppantPorosity, ProppantPermeability > >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    if( auto * const ptr = dynamic_cast< ProppantSolid< ProppantPorosity, ProppantPermeability > * >( &constitutiveRelation ) )
    {
      lambda( *ptr );
    }
    else
    {
      GEOS_ERROR( "ConstitutivePassThru< ProppantSolid >::execute failed on constitutive relation "
                  << constitutiveRelation.getDataContext() << " with type "
                  << LvArray::system::demangleType( constitutiveRelation ) );
    }
  }
};


/**
 * Specialization for all CoupledSolid models.
 */
template<>
struct ConstitutivePassThru< CoupledSolidBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< CompressibleSolid< PressurePorosity, ConstantPermeability >,
                                 CompressibleSolid< PressurePorosity, CarmanKozenyPermeability >,
                                 CompressibleSolid< PressurePorosity, ExponentialDecayPermeability >,
                                 CompressibleSolid< PressurePorosity, ParallelPlatesPermeability >,
                                 CompressibleSolid< PressurePorosity, PressurePermeability >,
                                 CompressibleSolid< PressurePorosity, SlipDependentPermeability >,
                                 CompressibleSolid< PressurePorosity, WillisRichardsPermeability >,
                                 PorousSolid< DruckerPragerExtended >,
                                 PorousSolid< ModifiedCamClay >,
                                 PorousSolid< DelftEgg >,
                                 PorousSolid< DruckerPrager >,
                                 PorousSolid< DuvautLionsSolid< DruckerPrager > >,
                                 PorousSolid< DuvautLionsSolid< DruckerPragerExtended > >,
                                 PorousSolid< DuvautLionsSolid< ModifiedCamClay > >,
                                 PorousSolid< ElasticIsotropic >,
                                 PorousSolid< ElasticTransverseIsotropic >,
                                 PorousSolid< ElasticIsotropicPressureDependent >,
                                 PorousSolid< ElasticOrthotropic >,
                                 PorousSolid< DamageSpectral< ElasticIsotropic > >,
                                 PorousSolid< DamageVolDev< ElasticIsotropic > >,
                                 PorousSolid< Damage< ElasticIsotropic > > >::execute( constitutiveRelation,
                                                                                       std::forward< LAMBDA >( lambda ) );
  }

  template< typename LAMBDA >
  static void execute( ConstitutiveBase const & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< CompressibleSolid< PressurePorosity, ConstantPermeability >,
                                 CompressibleSolid< PressurePorosity, CarmanKozenyPermeability >,
                                 CompressibleSolid< PressurePorosity, ExponentialDecayPermeability >,
                                 CompressibleSolid< PressurePorosity, ParallelPlatesPermeability >,
                                 CompressibleSolid< PressurePorosity, PressurePermeability >,
                                 CompressibleSolid< PressurePorosity, SlipDependentPermeability >,
                                 CompressibleSolid< PressurePorosity, WillisRichardsPermeability >,
                                 PorousSolid< DruckerPragerExtended >,
                                 PorousSolid< ModifiedCamClay >,
                                 PorousSolid< DelftEgg >,
                                 PorousSolid< DruckerPrager >,
                                 PorousSolid< DuvautLionsSolid< DruckerPrager > >,
                                 PorousSolid< DuvautLionsSolid< DruckerPragerExtended > >,
                                 PorousSolid< DuvautLionsSolid< ModifiedCamClay > >,
                                 PorousSolid< ElasticIsotropic >,
                                 PorousSolid< ElasticTransverseIsotropic >,
                                 PorousSolid< ElasticIsotropicPressureDependent >,
                                 PorousSolid< ElasticOrthotropic >,
                                 PorousSolid< DamageSpectral< ElasticIsotropic > >,
                                 PorousSolid< DamageVolDev< ElasticIsotropic > >,
                                 PorousSolid< Damage< ElasticIsotropic > > >::execute( constitutiveRelation,
                                                                                       std::forward< LAMBDA >( lambda ) );
  }
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_ */
