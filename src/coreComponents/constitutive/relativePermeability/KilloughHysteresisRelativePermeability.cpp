//
// Created by root on 11/9/22.
//

#include "KilloughHysteresisRelativePermeability.hpp"


namespace geosx
{

namespace constitutive
{

KilloughHysteresisRelativePermeability::KilloughHysteresisRelativePermeabilityKernel::KilloughHysteresisRelativePermeabilityKernel( const geosx::real64 & jerauldParam_a,
                                                                                                                                    const geosx::real64 & jerauldParam_b,
                                                                                                                                    const geosx::real64 & killoughCurvatureParam,
                                                                                                                                    const arrayView1d< const geosx::real64 > & landParam,
                                                                                                                                    const arrayView1d< const geosx::real64 > & drainageMinPhaseVolFraction,
                                                                                                                                    const arrayView1d< const geosx::real64 > & imbibitionMinPhaseVolFraction,
                                                                                                                                    const arrayView1d< const geosx::real64 > & drainageRelPermEndPoint,
                                                                                                                                    const arrayView1d< const geosx::real64 > & imbibitionRelPermEndPoint,
                                                                                                                                    const arrayView1d< const geosx::real64 > & drainageMaxPhaseVolFraction,
                                                                                                                                    const arrayView1d< const geosx::real64 > & imbibitionMaxPhaseVolFraction,
                                                                                                                                    const arrayView3d< geosx::real64,
                                                                                                                                                       relperm::USD_RELPERM > & phaseTrapppedVolFrac )
  :
  KilloughHysteresis::KernelKilloughHysteresisBase( jerauldParam_a,
                                                    jerauldParam_b,
                                                    landParam,
                                                    drainageMinPhaseVolFraction,
                                                    imbibitionMinPhaseVolFraction,
                                                    drainageRelPermEndPoint,
                                                    imbibitionRelPermEndPoint,
                                                    drainageMaxPhaseVolFraction,
                                                    imbibitionMaxPhaseVolFraction,
                                                    phaseTrapppedVolFrac
                                                    ),
  m_killoughCurvatureParam( killoughCurvatureParam )
{ };


KilloughHysteresisRelativePermeability::KilloughHysteresisRelativePermeability( const std::string & name,
                                                                                geosx::dataRepository::Group * const parent )
  : KilloughHysteresis( name, parent )
{

  void KilloughHysteresisRelativePermeability::postProcessInput()
  {

    KilloughHysteresis::postProcessInput();


  }


//TODO
//REGISTER_CATALOG_ENTRY()

}//end namespace constitutive
}//end namespace geosx
