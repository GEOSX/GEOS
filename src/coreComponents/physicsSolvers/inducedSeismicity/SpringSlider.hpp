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



#ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_SPRINGSLIDER_HPP
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_SPRINGSLIDER_HPP

#include "physicsSolvers/inducedSeismicity/ImplicitQDRateAndState.hpp"

namespace geos
{

template< typename RSSOLVER_TYPE = ImplicitQDRateAndState >
class SpringSlider : public RSSOLVER_TYPE
{
public:

  SpringSlider() = delete;

  SpringSlider( const string & name,
                dataRepository::Group * const parent );

  /// Destructor
  virtual ~SpringSlider() override;

  static string catalogName() { return RSSOLVER_TYPE::derivedSolverPrefix() + "SpringSlider"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override;

  struct viewKeyStruct : public RSSOLVER_TYPE::viewKeyStruct
  {};

private:

  virtual real64 updateStresses( real64 const & time_n,
                                 real64 const & dt,
                                 const int cycleNumber,
                                 DomainPartition & domain ) const override final;

  class SpringSliderParameters
  {
public:

    GEOS_HOST_DEVICE
    SpringSliderParameters( real64 const normalTraction, real64 const a, real64 const b, real64 const Dc ):
      tauRate( 1e-4 ),
      springStiffness( 0.0 )
    {
      real64 const criticalStiffness = normalTraction * (b - a) / Dc;
      springStiffness = 0.9 * criticalStiffness;
    }

    /// Default copy constructor
    SpringSliderParameters( SpringSliderParameters const & ) = default;

    /// Default move constructor
    SpringSliderParameters( SpringSliderParameters && ) = default;

    /// Deleted default constructor
    SpringSliderParameters() = delete;

    /// Deleted copy assignment operator
    SpringSliderParameters & operator=( SpringSliderParameters const & ) = delete;

    /// Deleted move assignment operator
    SpringSliderParameters & operator=( SpringSliderParameters && ) =  delete;

    real64 tauRate;

    real64 springStiffness;
  };
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_SPRINGSLIDER_HPP */
