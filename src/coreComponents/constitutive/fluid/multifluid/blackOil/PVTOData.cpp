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

#include "PVTOData.hpp"

namespace geos
{

namespace constitutive
{

PVTOData::KernelWrapper PVTOData::createKernelWrapper() const
{
  return PVTOData::KernelWrapper( Rs.toViewConst(),
                                  bubblePressure.toViewConst(),
                                  saturatedBo.toViewConst(),
                                  saturatedViscosity.toViewConst(),
                                  undersaturatedPressure2d.toViewConst(),
                                  undersaturatedBo2d.toViewConst(),
                                  undersaturatedViscosity2d.toViewConst(),
                                  surfaceMassDensity.toViewConst(),
                                  surfaceMoleDensity.toViewConst());
}

PVTOData::KernelWrapper::KernelWrapper( arrayView1d< real64 const > const & Rs,
                                        arrayView1d< real64 const > const & bubblePressure,
                                        arrayView1d< real64 const > const & saturatedBo,
                                        arrayView1d< real64 const > const & saturatedViscosity,
                                        arrayView2d< real64 const > const & undersaturatedPressure,
                                        arrayView2d< real64 const > const & undersaturatedBo,
                                        arrayView2d< real64 const > const & undersaturatedViscosity,
                                        arrayView1d< real64 const > const & surfaceMassDensity,
                                        arrayView1d< real64 const > const & surfaceMoleDensity )
  :
  m_Rs( Rs ),
  m_bubblePressure( bubblePressure ),
  m_saturatedBo( saturatedBo ),
  m_saturatedViscosity( saturatedViscosity ),
  m_undersaturatedPressure2d( undersaturatedPressure ),
  m_undersaturatedBo2d( undersaturatedBo ),
  m_undersaturatedViscosity2d( undersaturatedViscosity ),
  m_surfaceMassDensity( surfaceMassDensity ),
  m_surfaceMoleDensity( surfaceMoleDensity )
{}

} //namespace constitutive

} //namespace geos
