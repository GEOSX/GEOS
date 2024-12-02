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

/**
 * @file TableRelativePermeabilityHysteresis.hpp
 */

#ifndef GEOS_KILLOUGHHYSTERESIS_HPP
#define GEOS_KILLOUGHHYSTERESIS_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "functions/TableFunction.hpp"

#include "relativePermeability/Layouts.hpp"
#include "capillaryPressure/Layouts.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

/***
 * @brief KilloughHysteresis is designed to hold Killough hystereis model parameters and
 *        be in charge of all compuration related to this model (trapped Saturation,Land Coefficient?)
 */



//should be up to constitutiveBase or some new SCALConstitutiveBase but for now let's POC on relativePermeabilityBase
class KilloughHysteresis
{
public:

  static constexpr real64 minScriMinusScrd = 1e-12;
  /// To avoid frequent changes from drainage to imbibition and vice versa, we use this buffer
  static constexpr real64 flowReversalBuffer = 1e-12;

  struct PhaseWettability
  {
    enum : integer
    {
      WETTING = 0,
      NONWETTING = 1
    };
  };

  /**
   * @brief  struct to represent hysteresis curves (relperm or capillary pressure)
   * whether they are wetting or non wetting, storing key point as pairs of
   * saturations and value, being either the relperm value (S,kr) or capillary pressure value (S,pc).
   * this way we can tell wetting curve from non wetting by the ordering of drainage/imbibition key values.
   * Indeed if imbibition comes at lower saturation than drainage then it is wetting curve, on the opposite
   * if drainage comes before imbibition this is a non-wetting hysteresis. This is completed by an opposite
   * point that is either the connate wetting saturation Swc or the maximum non wetting saturation Sgmax.
   * @param  oppositeBoundPhaseVolFraction represents either Swc or Sgmax depending if a wetting curve or nonwetting is described
   * @param imbibitionExtremaPhaseVolFraction represents in wetting case the imibibition max and in non-wetting the imbibition residual
   * @param drainageExtremaPhaseVolFraction represents in wetting case the drainage max and in non-wetting the drainage residual
   * @param oppositeBoundSCALValue represents the associate relperm or capillary pressure value
   * @param imbibitionExtremaSCALValue represents the associate relperm or capillary pressure value
   * @param drainageExtremaSCALValue represents the associate relperm or capillary pressure value
   */

  struct HysteresisCurve
  {
    real64 oppositeBoundPhaseVolFraction = -1.;
    real64 imbibitionExtremaPhaseVolFraction = -1.;
    real64 drainageExtremaPhaseVolFraction = -1.;

    real64 oppositeBoundSCALValue = -1.;
    real64 imbibitionExtremaSCALValue = -1.;
    real64 drainageExtremaSCALValue = -1.;

    HysteresisCurve() = default;

    HysteresisCurve( std::pair< real64, real64 > const & opp, std::pair< real64, real64 > const & imbE, std::pair< real64, real64 > const & drainE )
    {
      setPoints( opp, imbE, drainE );
    }

    void setPoints( std::pair< real64, real64 > const & opp, std::pair< real64, real64 > const & imbE, std::pair< real64, real64 > const & drainE )
    {
      oppositeBoundPhaseVolFraction = opp.first;
      imbibitionExtremaPhaseVolFraction = imbE.first;
      drainageExtremaPhaseVolFraction = drainE.first;

      oppositeBoundSCALValue = opp.second;
      imbibitionExtremaSCALValue = imbE.second;
      drainageExtremaSCALValue = drainE.second;
    }

    //tODO (jacques) check if relevant to invert relation with same SCAL value // might be misleading for kr
    HysteresisCurve toWetting() const
    {
        if(!isWetting())
            return HysteresisCurve({1.-oppositeBoundPhaseVolFraction,oppositeBoundSCALValue},
                                   {1.- imbibitionExtremaPhaseVolFraction,imbibitionExtremaSCALValue},
                                   {1.-drainageExtremaPhaseVolFraction,drainageExtremaSCALValue});
        else
            return *this;
    }

      HysteresisCurve toNonWetting() const
      {
          if(isWetting())
              return HysteresisCurve({1.-oppositeBoundPhaseVolFraction,oppositeBoundSCALValue},
                                     {1.-imbibitionExtremaPhaseVolFraction,imbibitionExtremaSCALValue},
                                     {1.-drainageExtremaPhaseVolFraction,drainageExtremaSCALValue});
          else
              return *this;
      }

    bool isWetting() const
    {
      return ((drainageExtremaPhaseVolFraction > oppositeBoundPhaseVolFraction) ? PhaseWettability::WETTING : PhaseWettability::NONWETTING) == PhaseWettability::WETTING;
    }
    bool isZero() const
    {
      return (oppositeBoundPhaseVolFraction <= 0.0) && (imbibitionExtremaPhaseVolFraction <= 0.0) && (drainageExtremaPhaseVolFraction <= 0.0);
    }

  };

//  void setRelPermParameters( real64 const & jerauldA, real64 const & jerauldB, real64 const & relpermCurv );

  static std::string catalogName() { return "KilloughHysteresis"; }

  static void postProcessInput(real64 const &jerauldParam_a, real64 const &jerauldParam_b,
                               real64 const &killoughCurvatureParamRelPerm,
                               real64 const &killoughCurvatureParamPc);

  GEOS_HOST_DEVICE
  static void computeLandCoefficient( HysteresisCurve const & hcruve, real64 & landParam );
  /**
   * @brief Function computing the trapped critical phase volume fraction
   * @param[in] hcurve the hysteresis curve to be used and dispatched on
   * @param[in] Shy the max historical phase volume fraction
   * @param[in] landParam Land trapping parameter
   * @param[in] jerauldParam_a jerauld expononent
   * @param[in] jerauldParam_b jerauld expononent
   * @param[out] Scrt the trapped critical phase volume fraction
   */
  GEOS_HOST_DEVICE
  static void computeTrappedCriticalPhaseVolFraction( HysteresisCurve const & hcurve,
                                                      real64 const & Shy,
                                                      real64 const & landParam,
                                                      real64 const & jerauldParam_a,
                                                      real64 const & jerauldParam_b,
                                                      real64 & Scrt );

  struct viewKeyStruct
  {
    static constexpr char const * jerauldParameterAString() { return "jerauldParameterA"; }
    static constexpr char const * jerauldParameterBString() { return "jerauldParameterB"; }

    static constexpr char const * killoughCurvatureParameterRelPermString() { return "killoughCurvatureParameterRelPerm"; }
      static constexpr char const * killoughCurvatureParameterPcString() { return "killoughCurvatureParameterPc"; }
  };

};

}

}

#endif //GEOS_KILLOUGHHYSTERESIS_HPP
