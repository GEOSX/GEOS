/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/***
 *  @file KilloughHysteresis.cpp
 */

#include "KilloughHysteresis.hpp"


namespace geos
{

using namespace dataRepository;

namespace constitutive
{


void KilloughHysteresis::postProcessInput(real64 const &jerauldParam_a, real64 const &jerauldParam_b,
                                          real64 const &killoughCurvatureParamRelPerm,
                                          real64 const &killoughCurvatureParamPc)
{
  GEOS_THROW_IF( jerauldParam_a < 0,
                  GEOS_FMT( "{}: the parameter {} must be positive",
                             catalogName(),
                             viewKeyStruct::jerauldParameterAString() ),
                  InputError );

  GEOS_THROW_IF( jerauldParam_b < 0,
                  GEOS_FMT( "{}: the paramater {} must be postitive",
                             catalogName(),
                             viewKeyStruct::jerauldParameterBString() ),
                  InputError );

  GEOS_THROW_IF( killoughCurvatureParamRelPerm < 0,
                  GEOS_FMT( "{}: the paramater {} must be postitive",
                             catalogName(),
                             viewKeyStruct::killoughCurvatureParameterRelPermString() ),
                  InputError );

  GEOS_THROW_IF( killoughCurvatureParamPc < 0,
                    GEOS_FMT( "{}: the paramater {} must be postitive",
                               catalogName(),
                               viewKeyStruct::killoughCurvatureParameterPcString() ),
                    InputError );

}



//TODO
void KilloughHysteresis::computeLandCoefficient( KilloughHysteresis::HysteresisCurve const & hcurve,
                                                 real64 & landParam )
{

  // Note: for simplicity, the notations are taken from IX documentation (although this breaks our phaseVolFrac naming convention)

  // Step 1: Land parameter for the wetting phase
  if( hcurve.isWetting() )
  {
    real64 const Scrd = hcurve.oppositeBoundPhaseVolFraction;
    real64 const Smxd = hcurve.drainageExtremaPhaseVolFraction;
    real64 const Smxi = hcurve.imbibitionExtremaPhaseVolFraction;
    real64 const Swc = Scrd;
    GEOS_THROW_IF(  (Smxi - Smxd) > 0,
                     GEOS_FMT( "{}: For wetting phase hysteresis, imbibition end-point saturation Smxi( {} ) must be smaller than the drainage saturation end-point Smxd( {} ).\n"
                                "Crossing relative permeability curves.\n",
                                catalogName(),
                                Smxi,
                                Smxd ),
                     InputError );

    landParam = ( Smxd - Swc ) / LvArray::math::max( KilloughHysteresis::minScriMinusScrd, ( Smxd - Smxi ) ) - 1.0;
  }
  else
  // Step 2: Land parameter for the non-wetting phase

  {
    real64 const Smx =  hcurve.oppositeBoundPhaseVolFraction;
    real64 const Scrd = hcurve.drainageExtremaPhaseVolFraction;
    real64 const Scri = hcurve.imbibitionExtremaPhaseVolFraction;
    GEOS_THROW_IF( (Scrd - Scri) > 0,
                    GEOS_FMT( "{}: For non-wetting phase hysteresis, drainage trapped saturation Scrd( {} ) must be smaller than the imbibition saturation Scri( {} ).\n"
                               "Crossing relative permeability curves.\n",
                               catalogName(),
                               Scrd,
                               Scri ),
                    InputError );

    landParam = ( Smx - Scrd ) / LvArray::math::max( KilloughHysteresis::minScriMinusScrd, ( Scri - Scrd ) ) - 1.0;
  }
}

GEOS_HOST_DEVICE
void
KilloughHysteresis::
  computeTrappedCriticalPhaseVolFraction( HysteresisCurve const & hcurve,
                                          real64 const & Shy,
                                          real64 const & landParam,
                                          real64 const & jerauldParam_a,
                                          real64 const & jerauldParam_b,
                                          real64 & Scrt )
{

  if( hcurve.isWetting())
  {
    //unpack values
    real64 const Smxd = hcurve.drainageExtremaPhaseVolFraction;
    real64 const Swc = hcurve.oppositeBoundPhaseVolFraction;

    real64 const A = 1 + jerauldParam_a * (Shy - Swc);
    real64 const numerator = Shy - Smxd;
    real64 const denom = A + landParam * pow((Smxd - Shy) / (Smxd - Swc), 1 + jerauldParam_b / landParam );
    Scrt = Smxd + numerator / denom;
  }
  else
  {
    //unpack values
    real64 const Scrd = hcurve.drainageExtremaPhaseVolFraction;
    real64 const Smx = hcurve.oppositeBoundPhaseVolFraction;

    real64 const A = 1 + jerauldParam_a * (Smx - Shy);
    real64 const numerator = Shy - Scrd;
    real64 const denom = A + landParam * pow((Shy - Scrd) / (Smx - Scrd), 1 + jerauldParam_b / landParam );
    Scrt = LvArray::math::max( 0.0,
                               Scrd + numerator / denom );           // trapped critical saturation from equation 2.162
  }

}

}//end namespace
}//end namespace
