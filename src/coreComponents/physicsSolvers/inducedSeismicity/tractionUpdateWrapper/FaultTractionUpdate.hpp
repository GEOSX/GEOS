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
 * @file FaultTractionUpdate.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_FAULTTRACTIONUPDATE_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_FAULTTRACTIONUPDATE_HPP


#include "physicsSolvers/inducedSeismicity/tractionUpdateWrapper/FaultTractionUpdateBase.hpp"


namespace geos
{

namespace inducedSeismicity
{

template< typename ... SOLVER_TYPE >
class FaultTractionUpdate : public FaultTractionUpdateBase
{
public:

  FaultTractionUpdate( SOLVER_TYPE * ... solvers ):
    m_solvers( std::make_tuple( solvers ... ) )
  {}

  virtual ~FaultTractionUpdate() = default;

protected:

  std::tuple< SOLVER_TYPE * ... > m_solvers;

};

}

}

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_FAULTTRACTIONUPDATE_HPP */
