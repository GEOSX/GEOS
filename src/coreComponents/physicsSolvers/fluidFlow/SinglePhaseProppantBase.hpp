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
 * @file SinglePhaseProppantBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEPROPPANTBASE_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEPROPPANTBASE_HPP_

#include "SinglePhaseBase.hpp"

namespace geos
{

class SinglePhaseProppantBase : public SinglePhaseBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseProppantBase( const string & name,
                           Group * const parent );

  SinglePhaseProppantBase() = delete;

  /// deleted copy constructor
  SinglePhaseProppantBase( SinglePhaseProppantBase const & ) = delete;

  /// default move constructor
  SinglePhaseProppantBase( SinglePhaseProppantBase && ) = default;

  /// deleted assignment operator
  SinglePhaseProppantBase & operator=( SinglePhaseProppantBase const & ) = delete;

  /// deleted move operator
  SinglePhaseProppantBase & operator=( SinglePhaseProppantBase && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseProppantBase();

  virtual void updateFluidModel( ObjectManagerBase & dataGroup ) const override;

  virtual void updatePorosityAndPermeability( SurfaceElementSubRegion & subRegion ) const override;

protected:

  virtual void validateConstitutiveModels( DomainPartition & domain ) const override;

  virtual FluidPropViews getFluidProperties( constitutive::ConstitutiveBase const & fluid ) const override;

private:
  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;


};
}
#endif /* GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEPROPPANTBASE_HPP_ */
