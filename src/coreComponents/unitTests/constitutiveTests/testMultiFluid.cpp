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

// Source includes
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "mainInterface/initialization.hpp"
#include "functions/FunctionManager.hpp"
#include "mainInterface/GeosxState.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;
using namespace geos::constitutive::multifluid;
using namespace geos::dataRepository;
using namespace geos::constitutive::PVTProps;

/// Black-oil tables written into temporary files during testing

static const char * pvtoTableContent = "# Rs[sm3/sm3]\tPbub[Pa]\tBo[m3/sm3]\tVisc(Pa.s)\n"
                                       "\n"
                                       "  2\t            2000000\t    1.02\t    0.000975\n"
                                       "  5\t            5000000\t    1.03\t    0.00091\n"
                                       " 10\t            10000000\t1.04\t    0.00083\n"
                                       " 15\t            20000000\t1.05\t    0.000695\n"
                                       "                90000000\t1.03\t    0.000985  -- some line comment\n"
                                       " 30\t            30000000\t1.07\t    0.000594\n"
                                       " 40\t            40000000\t1.08\t    0.00051\n"
                                       "                50000000\t1.07\t    0.000549  -- another one\n"
                                       "                90000000\t1.06\t    0.00074\n"
                                       " 50\t            50000000.7\t1.09\t    0.000449\n"
                                       "                90000000.7\t1.08\t    0.000605";

static const char * pvtwTableContent = "#\tPref[Pa]\tBw[m3/sm3]\tCp[1/Pa]\t    Visc[Pa.s]\n"
                                       "\t30600000.1\t1.03\t\t0.00000000041\t0.0003";

/// Dead-oil tables written into temporary files during testing

static const char * pvdgTableContent = "# Pg(Pa) Bg(m3/sm3) Visc(Pa.s)\n"
                                       "3000000  0.04234  0.00001344\n"
                                       "6000000  0.02046  0.0000142\n"
                                       "9000000  0.01328  0.00001526\n"
                                       "12000000 0.00977  0.0000166\n"
                                       "15000000 0.00773  0.00001818\n"
                                       "18000000 0.006426 0.00001994\n"
                                       "21000000 0.005541 0.00002181\n"
                                       "24000000 0.004919 0.0000237\n"
                                       "27000000 0.004471 0.00002559\n"
                                       "29500000 0.004194 0.00002714\n"
                                       "31000000 0.004031 0.00002806\n"
                                       "33000000 0.00391  0.00002832\n"
                                       "53000000 0.003868 0.00002935";

static const char * pvdoTableContent = "#P[Pa] Bo[m3/sm3] Visc(Pa.s)\n"
                                       "10000000.0 1.23331 0.00015674\n"
                                       "12500000.0 1.21987 0.00016570\n"
                                       "15000000.0 1.20802 0.00017445\n"
                                       "20000000.0 1.18791 0.00019143\n"
                                       "25000000.0 1.17137 0.00020779\n"
                                       "30000000.0 1.15742 0.00022361\n"
                                       "33200000.3 1.14946 0.00023359\n"
                                       "35000000.0 1.14543 0.00023894\n"
                                       "40000000.0 1.13498 0.00025383\n"
                                       "50000000.0 1.11753 0.00028237\n"
                                       "60000000.0 1.10346 0.00030941\n"
                                       "70000000.0 1.09180 0.00033506\n"
                                       "80000000.0 1.08194 0.00035945\n"
                                       "90000000.0 1.07347 0.00038266\n"
                                       "95000000.0 1.06966 0.00039384\n"
                                       "100000000.0 1.06610 0.00040476\n"
                                       "110000000.0 1.05961 0.00042584\n"
                                       "112500000.0 1.05811 0.00043096\n"
                                       "115000000.0 1.05665 0.00043602\n"
                                       "117500000.0 1.05523 0.00044102\n"
                                       "120000000.0 1.05385 0.00044596\n";

static const char * pvdwTableContent = "# Pref[Pa] Bw[m3/sm3] Cp[1/Pa]     Visc[Pa.s]\n"
                                       " 30600000.1 1.03  0.00000000041 0.0003";

// CO2-brine model

static const char * pvtLiquidPhillipsTableContent = "DensityFun PhillipsBrineDensity 1e6 1.5e7 5e4 367.15 369.15 1 0.2\n"
                                                    "ViscosityFun PhillipsBrineViscosity 0.1";

// the last are set relatively high (1e-4) to increase derivative value and check it properly
static const char * pvtLiquidEzrokhiTableContent = "DensityFun EzrokhiBrineDensity 2.01e-6 -6.34e-7 1e-4\n"
                                                   "ViscosityFun EzrokhiBrineViscosity 2.42e-7 0 1e-4";

static const char * pvtGasTableContent = "DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                         "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 367.15 369.15 1";

static const char * co2FlashTableContent = "FlashModel CO2Solubility 1e6 1.5e7 5e4 367.15 369.15 1 0.15";

void testNumericalDerivatives( MultiFluidBase & fluid,
                               Group & parent,
                               real64 const P,
                               real64 const T,
                               arraySlice1d< real64 > const & compositionInput,
                               real64 const perturbParameter,
                               bool usePVTPackage,
                               real64 const relTol,
                               real64 const absTol = std::numeric_limits< real64 >::max() )
{
  using Deriv = multifluid::DerivativeOffset;

  integer const NC = fluid.numFluidComponents();
  integer const NP = fluid.numFluidPhases();
  integer const NDOF = NC+2;

  // Copy input values into an array with expected layout
  array2d< real64, compflow::LAYOUT_COMP > compositionValues( 1, NC );
  for( integer i = 0; i < NC; ++i )
  {
    compositionValues[0][i] = compositionInput[i];
  }
  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition = compositionValues[0];

  auto const & components = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  auto const & phases     = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );

  // create a clone of the fluid to run updates on
  std::unique_ptr< ConstitutiveBase > fluidCopyPtr = fluid.deliverClone( "fluidCopy", &parent );
  MultiFluidBase & fluidCopy = dynamicCast< MultiFluidBase & >( *fluidCopyPtr );

  fluid.allocateConstitutiveData( fluid.getParent(), 1 );
  fluidCopy.allocateConstitutiveData( fluid.getParent(), 1 );

  // extract data views from both fluids
  #define GET_FLUID_DATA( FLUID, TRAIT ) \
    FLUID.getReference< TRAIT::type >( TRAIT::key() )[0][0]

  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 > phaseFrac {
    GET_FLUID_DATA( fluid, fields::multifluid::phaseFraction ),
    GET_FLUID_DATA( fluid, fields::multifluid::dPhaseFraction )
  };

  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 > phaseDens {
    GET_FLUID_DATA( fluid, fields::multifluid::phaseDensity ),
    GET_FLUID_DATA( fluid, fields::multifluid::dPhaseDensity )
  };

  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 > phaseVisc {
    GET_FLUID_DATA( fluid, fields::multifluid::phaseViscosity ),
    GET_FLUID_DATA( fluid, fields::multifluid::dPhaseViscosity )
  };

  MultiFluidVarSlice< real64, 2, USD_PHASE_COMP - 2, USD_PHASE_COMP_DC - 2 > phaseCompFrac {
    GET_FLUID_DATA( fluid, fields::multifluid::phaseCompFraction ),
    GET_FLUID_DATA( fluid, fields::multifluid::dPhaseCompFraction )
  };

  MultiFluidVarSlice< real64, 0, USD_FLUID - 2, USD_FLUID_DC - 2 > totalDens {
    GET_FLUID_DATA( fluid, fields::multifluid::totalDensity ),
    GET_FLUID_DATA( fluid, fields::multifluid::dTotalDensity )
  };

  auto const & phaseFracCopy     = GET_FLUID_DATA( fluidCopy, fields::multifluid::phaseFraction );
  auto const & phaseDensCopy     = GET_FLUID_DATA( fluidCopy, fields::multifluid::phaseDensity );
  auto const & phaseViscCopy     = GET_FLUID_DATA( fluidCopy, fields::multifluid::phaseViscosity );
  auto const & phaseCompFracCopy = GET_FLUID_DATA( fluidCopy, fields::multifluid::phaseCompFraction );
  auto totalDensCopy = [&]() { return GET_FLUID_DATA( fluidCopy, fields::multifluid::totalDensity ); };

#undef GET_FLUID_DATA

  // set the original fluid state to current
  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    fluidWrapper.update( 0, 0, P, T, composition );
  } );

  // now perturb variables and update the copied fluid's state
  constitutive::constitutiveUpdatePassThru( fluidCopy, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    // to be able to use the checkDerivative utility function, we have to invert the layout
    auto dPhaseFrac     = invertLayout( phaseFrac.derivs.toSliceConst(), NP, NDOF );
    auto dPhaseDens     = invertLayout( phaseDens.derivs.toSliceConst(), NP, NDOF );
    auto dPhaseVisc     = invertLayout( phaseVisc.derivs.toSliceConst(), NP, NDOF );
    auto dTotalDens     = invertLayout( totalDens.derivs.toSliceConst(), NDOF );
    auto dPhaseCompFrac = invertLayout( phaseCompFrac.derivs.toSliceConst(), NP, NC, NDOF );

    // update pressure and check derivatives
    {
      real64 const dP = perturbParameter * (P + perturbParameter);
      fluidWrapper.update( 0, 0, P + dP, T, composition );

      checkDerivative( phaseFracCopy.toSliceConst(), phaseFrac.value.toSliceConst(), dPhaseFrac[Deriv::dP].toSliceConst(),
                       dP, relTol, absTol, "phaseFrac", "Pres", phases );
      checkDerivative( phaseDensCopy.toSliceConst(), phaseDens.value.toSliceConst(), dPhaseDens[Deriv::dP].toSliceConst(),
                       dP, relTol, absTol, "phaseDens", "Pres", phases );
      checkDerivative( phaseViscCopy.toSliceConst(), phaseVisc.value.toSliceConst(), dPhaseVisc[Deriv::dP].toSliceConst(),
                       dP, relTol, absTol, "phaseVisc", "Pres", phases );
      checkDerivative( totalDensCopy(), totalDens.value, dTotalDens[Deriv::dP],
                       dP, relTol, absTol, "totalDens", "Pres" );
      checkDerivative( phaseCompFracCopy.toSliceConst(), phaseCompFrac.value.toSliceConst(), dPhaseCompFrac[Deriv::dP].toSliceConst(),
                       dP, relTol, absTol, "phaseCompFrac", "Pres", phases, components );
    }

    // update temperature and check derivatives
    {
      real64 const dT = perturbParameter * (T + perturbParameter);
      fluidWrapper.update( 0, 0, P, T + dT, composition );

      checkDerivative( phaseFracCopy.toSliceConst(), phaseFrac.value.toSliceConst(), dPhaseFrac[Deriv::dT].toSliceConst(),
                       dT, relTol, absTol, "phaseFrac", "Temp", phases );
      checkDerivative( phaseDensCopy.toSliceConst(), phaseDens.value.toSliceConst(), dPhaseDens[Deriv::dT].toSliceConst(),
                       dT, relTol, absTol, "phaseDens", "Temp", phases );
      checkDerivative( phaseViscCopy.toSliceConst(), phaseVisc.value.toSliceConst(), dPhaseVisc[Deriv::dT].toSliceConst(),
                       dT, relTol, absTol, "phaseVisc", "Temp", phases );
      checkDerivative( totalDensCopy(), totalDens.value, dTotalDens[Deriv::dT],
                       dT, relTol, absTol, "totalDens", "Temp" );
      checkDerivative( phaseCompFracCopy.toSliceConst(), phaseCompFrac.value.toSliceConst(), dPhaseCompFrac[Deriv::dT].toSliceConst(),
                       dT, relTol, absTol, "phaseCompFrac", "Temp", phases, components );
    }

    array2d< real64, compflow::LAYOUT_COMP > compNew( 1, NC );
    for( integer jc = 0; jc < NC; ++jc )
    {
      real64 const dC = perturbParameter * ( composition[jc] + perturbParameter );
      for( integer ic = 0; ic < NC; ++ic )
      {
        compNew[0][ic] = composition[ic];
      }
      compNew[0][jc] += dC;

      // Note: in PVTPackage, derivatives are obtained with finite-difference approx **with normalization of the comp fraction**
      //       The component fraction is perturbed (just as above), and then all the component fractions are normalized (as below)
      //       But, in the native DO model and in CO2BrinePhillips, derivatives are computed analytically, which results in different
      //       derivatives wrt component fractions--although the derivatives wrt component densities obtained with the chain rule
      //       in the solver will be very similar (see discussion on PR #1325 on GitHub).
      //
      //       Since both approaches--FD approximation of derivatives with normalization, and analytical derivatives--are correct,
      //       we have to support both when we check the intermediate derivatives wrt component fractions below. Therefore, if the
      //       PVTPackage is used, then we normalize the perturbed component fractions before taking the FD approx. If the native
      //       DO or CO2-brine models are used, we skip the normalization below.
      if( usePVTPackage )
      {
        // renormalize
        real64 sum = 0.0;
        for( integer ic = 0; ic < NC; ++ic )
        {
          sum += compNew[0][ic];
        }
        for( integer ic = 0; ic < NC; ++ic )
        {
          compNew[0][ic] /= sum;
        }
      }

      fluidWrapper.update( 0, 0, P, T, compNew[0] );

      string const var = "compFrac[" + components[jc] + "]";
      checkDerivative( phaseFracCopy.toSliceConst(), phaseFrac.value.toSliceConst(), dPhaseFrac[Deriv::dC+jc].toSliceConst(),
                       dC, relTol, absTol, "phaseFrac", var, phases );
      checkDerivative( phaseDensCopy.toSliceConst(), phaseDens.value.toSliceConst(), dPhaseDens[Deriv::dC+jc].toSliceConst(),
                       dC, relTol, absTol, "phaseDens", var, phases );
      checkDerivative( phaseViscCopy.toSliceConst(), phaseVisc.value.toSliceConst(), dPhaseVisc[Deriv::dC+jc].toSliceConst(),
                       dC, relTol, absTol, "phaseVisc", var, phases );
      checkDerivative( totalDensCopy(), totalDens.value, dTotalDens[Deriv::dC+jc],
                       dC, relTol, absTol, "totalDens", var );
      checkDerivative( phaseCompFracCopy.toSliceConst(), phaseCompFrac.value.toSliceConst(), dPhaseCompFrac[Deriv::dC+jc].toSliceConst(),
                       dC, relTol, absTol, "phaseCompFrac", var, phases, components );
    }
  } );
}

void testValuesAgainstPreviousImplementation( CO2BrinePhillipsFluid::KernelWrapper const & wrapper,
                                              real64 const P,
                                              real64 const T,
                                              arraySlice1d< real64 const > const & compositionInput,
                                              real64 const & savedTotalDensity,
                                              real64 const & savedGasPhaseFrac,
                                              real64 const & savedWaterDens,
                                              real64 const & savedGasDens,
                                              real64 const & savedWaterMassDens,
                                              real64 const & savedGasMassDens,
                                              real64 const & savedWaterVisc,
                                              real64 const & savedGasVisc,
                                              real64 const & savedWaterPhaseGasComp,
                                              real64 const & savedWaterPhaseWaterComp,
                                              real64 const relTol )
{
  integer constexpr numPhase = 2;
  integer constexpr numComp  = 2;
  integer constexpr numDof   = numComp + 2;

  // Copy input values into an array with expected layout
  array2d< real64, compflow::LAYOUT_COMP > compositionValues( 1, numComp );
  for( integer i = 0; i < numComp; ++i )
  {
    compositionValues[0][i] = compositionInput[i];
  }
  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition = compositionValues[0];

  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseFraction( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseFraction( 1, 1, numPhase, numDof );
  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseDensity( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseDensity( 1, 1, numPhase, numDof );
  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseMassDensity( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseMassDensity( 1, 1, numPhase, numDof );
  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseViscosity( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseViscosity( 1, 1, numPhase, numDof );
  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseEnthalpy( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseEnthalpy( 1, 1, numPhase, numDof );
  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseInternalEnergy( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseInternalEnergy( 1, 1, numPhase, numDof );
  StackArray< real64, 4, numComp *numPhase, LAYOUT_PHASE_COMP > phaseCompFraction( 1, 1, numPhase, numComp );
  StackArray< real64, 5, numDof *numComp *numPhase, LAYOUT_PHASE_COMP_DC > dPhaseCompFraction( 1, 1, numPhase, numComp, numDof );
  StackArray< real64, 2, 1, LAYOUT_FLUID > totalDensity( 1, 1 );
  StackArray< real64, 3, numDof, LAYOUT_FLUID_DC >  dTotalDensity( 1, 1, numDof );

  wrapper.compute( P, T, composition,
  {
    phaseFraction[0][0],
    dPhaseFraction[0][0]
  },
  {
    phaseDensity[0][0],
    dPhaseDensity[0][0]
  },
  {
    phaseMassDensity[0][0],
    dPhaseMassDensity[0][0]
  },
  {
    phaseViscosity[0][0],
    dPhaseViscosity[0][0]
  },
  {
    phaseEnthalpy[0][0],
    dPhaseEnthalpy[0][0]
  },
  {
    phaseInternalEnergy[0][0],
    dPhaseInternalEnergy[0][0]
  },
  {
    phaseCompFraction[0][0],
    dPhaseCompFraction[0][0]
  },
  {
    totalDensity[0][0],
    dTotalDensity[0][0]
  } );

  checkRelativeError( totalDensity[0][0], savedTotalDensity, relTol );
  for( integer ip = 0; ip < numPhase; ++ip )
  {
    real64 const savedPhaseFrac = ( ip == 0 ) ? savedGasPhaseFrac : 1 - savedGasPhaseFrac;
    checkRelativeError( phaseFraction[0][0][ip], savedPhaseFrac, relTol );
    real64 const savedPhaseDens = ( ip == 0 ) ? savedGasDens : savedWaterDens;
    checkRelativeError( phaseDensity[0][0][ip], savedPhaseDens, relTol );
    real64 const savedPhaseMassDens = ( ip == 0 ) ? savedGasMassDens : savedWaterMassDens;
    checkRelativeError( phaseMassDensity[0][0][ip], savedPhaseMassDens, relTol );
    real64 const savedPhaseVisc = ( ip == 0 ) ? savedGasVisc : savedWaterVisc;
    checkRelativeError( phaseViscosity[0][0][ip], savedPhaseVisc, relTol );
    for( integer ic = 0; ic < numComp; ++ic )
    {
      real64 savedCompFrac = 0.0;
      if( ip == 0 )
      {
        savedCompFrac = ( ic == 0 ) ? 1 : 0;
      }
      else
      {
        savedCompFrac = ( ic == 0 ) ? savedWaterPhaseGasComp : savedWaterPhaseWaterComp;
      }
      checkRelativeError( phaseCompFraction[0][0][ip][ic], savedCompFrac, relTol );
    }
  }
}

MultiFluidBase & makeCompositionalFluid( string const & name, Group & parent )
{
  CompositionalMultiphaseFluidPVTPackage & fluid = parent.registerGroup< CompositionalMultiphaseFluidPVTPackage >( name );

  // TODO we should actually create a fake XML node with data, but this seemed easier...

  auto & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 4 );
  compNames[0] = "N2"; compNames[1] = "C10"; compNames[2] = "C20"; compNames[3] = "H20";

  auto & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 4 );
  molarWgt[0] = 28e-3; molarWgt[1] = 134e-3; molarWgt[2] = 275e-3; molarWgt[3] = 18e-3;

  auto & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  auto & eqnOfState = fluid.getReference< string_array >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::equationsOfStateString() );
  eqnOfState.resize( 2 );
  eqnOfState[0] = "PR"; eqnOfState[1] = "PR";

  auto & critPres = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::componentCriticalPressureString() );
  critPres.resize( 4 );
  critPres[0] = 34e5; critPres[1] = 25.3e5; critPres[2] = 14.6e5; critPres[3] = 220.5e5;

  auto & critTemp = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::componentCriticalTemperatureString() );
  critTemp.resize( 4 );
  critTemp[0] = 126.2; critTemp[1] = 622.0; critTemp[2] = 782.0; critTemp[3] = 647.0;

  auto & acFactor = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::componentAcentricFactorString() );
  acFactor.resize( 4 );
  acFactor[0] = 0.04; acFactor[1] = 0.443; acFactor[2] = 0.816; acFactor[3] = 0.344;

  fluid.postProcessInputRecursive();
  return fluid;
}

class CompositionalFluidTestBase : public ::testing::Test
{
public:
  CompositionalFluidTestBase():
    node(),
    parent( "parent", node )
  {}

protected:
  conduit::Node node;
  Group parent;
  MultiFluidBase * fluid;
};

class CompositionalFluidTest : public CompositionalFluidTestBase
{
public:
  CompositionalFluidTest()
  {
    parent.resize( 1 );
    fluid = &makeCompositionalFluid( "fluid", parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }
};

TEST_F( CompositionalFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d< real64 > comp( 4 );
  comp[0] = 0.099; comp[1] = 0.3; comp[2] = 0.6; comp[3] = 0.001;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;

  testNumericalDerivatives( *fluid, parent, P, T, comp, eps, true, relTol );
}

TEST_F( CompositionalFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d< real64 > comp( 4 );
  comp[0] = 0.099; comp[1] = 0.3; comp[2] = 0.6; comp[3] = 0.001;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-2;

  testNumericalDerivatives( *fluid, parent, P, T, comp, eps, true, relTol );
}

MultiFluidBase & makeLiveOilFluid( string const & name, Group * parent )
{
  BlackOilFluid & fluid = parent->registerGroup< BlackOilFluid >( name );

  string_array & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 3 );
  compNames[0] = "oil"; compNames[1] = "gas"; compNames[2] = "water";

  array1d< real64 > & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 3 );
  molarWgt[0] = 114e-3; molarWgt[1] = 16e-3; molarWgt[2] = 18e-3;

  string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  array1d< real64 > & surfaceDens = fluid.getReference< array1d< real64 > >( BlackOilFluidBase::viewKeyStruct::surfacePhaseMassDensitiesString() );
  surfaceDens.resize( 3 );
  surfaceDens[0] = 800.0; surfaceDens[1] = 0.9907; surfaceDens[2] = 1022.0;

  path_array & tableNames = fluid.getReference< path_array >( BlackOilFluidBase::viewKeyStruct::tableFilesString() );
  tableNames.resize( 3 );
  tableNames[0] = "pvto.txt"; tableNames[1] = "pvdg.txt"; tableNames[2] = "pvtw.txt";

  fluid.postProcessInputRecursive();
  return fluid;
}

MultiFluidBase & makeDeadOilFluid( string const & name, Group * parent )
{
  DeadOilFluid & fluid = parent->registerGroup< DeadOilFluid >( name );

  string_array & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 3 );
  compNames[0] = "oil"; compNames[1] = "water"; compNames[2] = "gas";

  array1d< real64 > & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 3 );
  molarWgt[0] = 114e-3; molarWgt[1] = 16e-3; molarWgt[2] = 18e-3;

  string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "water"; phaseNames[2] = "gas";

  array1d< real64 > & surfaceDens = fluid.getReference< array1d< real64 > >( BlackOilFluidBase::viewKeyStruct::surfacePhaseMassDensitiesString() );
  surfaceDens.resize( 3 );
  surfaceDens[0] = 800.0; surfaceDens[1] = 1022.0; surfaceDens[2] = 0.9907;

  path_array & tableNames = fluid.getReference< path_array >( BlackOilFluidBase::viewKeyStruct::tableFilesString() );
  tableNames.resize( 3 );
  tableNames[0] = "pvdo.txt"; tableNames[1] = "pvdw.txt"; tableNames[2] = "pvdg.txt";

  fluid.postProcessInputRecursive();
  return fluid;
}

MultiFluidBase & makeDeadOilFluidFromTable( string const & name, Group * parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables (PVDO, PVDG)

  // 1D table with linear interpolation
  integer const NaxisPVDO = 21;
  integer const NaxisPVDG = 13;

  array1d< real64_array > coordinatesPVDO;
  real64_array valuesPVDO_Bo( NaxisPVDO );
  real64_array valuesPVDO_visc( NaxisPVDO );
  coordinatesPVDO.resize( 1 );
  coordinatesPVDO[0].resize( NaxisPVDO );
  coordinatesPVDO[0][0] = 10000000.0; valuesPVDO_Bo[0] = 1.23331; valuesPVDO_visc[0] = 0.00015674;
  coordinatesPVDO[0][1] = 12500000.0; valuesPVDO_Bo[1] = 1.21987; valuesPVDO_visc[1] = 0.00016570;
  coordinatesPVDO[0][2] = 15000000.0; valuesPVDO_Bo[2] = 1.20802; valuesPVDO_visc[2] = 0.00017445;
  coordinatesPVDO[0][3] = 20000000.0; valuesPVDO_Bo[3] = 1.18791; valuesPVDO_visc[3] = 0.00019143;
  coordinatesPVDO[0][4] = 25000000.0; valuesPVDO_Bo[4] = 1.17137; valuesPVDO_visc[4] = 0.00020779;
  coordinatesPVDO[0][5] = 30000000.0; valuesPVDO_Bo[5] = 1.15742; valuesPVDO_visc[5] = 0.00022361;
  coordinatesPVDO[0][6] = 33200000.3; valuesPVDO_Bo[6] = 1.14946; valuesPVDO_visc[6] = 0.00023359;
  coordinatesPVDO[0][7] = 35000000.0; valuesPVDO_Bo[7] = 1.14543; valuesPVDO_visc[7] = 0.00023894;
  coordinatesPVDO[0][8] = 40000000.0; valuesPVDO_Bo[8] = 1.13498; valuesPVDO_visc[8] = 0.00025383;
  coordinatesPVDO[0][9] = 50000000.0; valuesPVDO_Bo[9] = 1.11753; valuesPVDO_visc[9] = 0.00028237;
  coordinatesPVDO[0][10] = 60000000.0; valuesPVDO_Bo[10] = 1.10346; valuesPVDO_visc[10] = 0.00030941;
  coordinatesPVDO[0][11] = 70000000.0; valuesPVDO_Bo[11] = 1.09180; valuesPVDO_visc[11] = 0.00033506;
  coordinatesPVDO[0][12] = 80000000.0; valuesPVDO_Bo[12] = 1.08194; valuesPVDO_visc[12] = 0.00035945;
  coordinatesPVDO[0][13] = 90000000.0; valuesPVDO_Bo[13] = 1.07347; valuesPVDO_visc[13] = 0.00038266;
  coordinatesPVDO[0][14] = 95000000.0; valuesPVDO_Bo[14] = 1.06966; valuesPVDO_visc[14] = 0.00039384;
  coordinatesPVDO[0][15] = 100000000.0; valuesPVDO_Bo[15] = 1.06610; valuesPVDO_visc[15] = 0.00040476;
  coordinatesPVDO[0][16] = 110000000.0; valuesPVDO_Bo[16] = 1.05961; valuesPVDO_visc[16] = 0.00042584;
  coordinatesPVDO[0][17] = 112500000.0; valuesPVDO_Bo[17] = 1.05811; valuesPVDO_visc[17] = 0.00043096;
  coordinatesPVDO[0][18] = 115000000.0; valuesPVDO_Bo[18] = 1.05665; valuesPVDO_visc[18] = 0.00043602;
  coordinatesPVDO[0][19] = 117500000.0; valuesPVDO_Bo[19] = 1.05523; valuesPVDO_visc[19] = 0.00044102;
  coordinatesPVDO[0][20] = 120000000.0; valuesPVDO_Bo[20] = 1.05385; valuesPVDO_visc[20] = 0.00044596;

  array1d< real64_array > coordinatesPVDG;
  real64_array valuesPVDG_Bg( NaxisPVDG );
  real64_array valuesPVDG_visc( NaxisPVDG );
  coordinatesPVDG.resize( 1 );
  coordinatesPVDG[0].resize( NaxisPVDG );
  coordinatesPVDG[0][0]  =  3000000; valuesPVDG_Bg[0]  = 0.04234;  valuesPVDG_visc[0] = 0.00001344;
  coordinatesPVDG[0][1]  =  6000000; valuesPVDG_Bg[1]  = 0.02046;  valuesPVDG_visc[1] = 0.0000142;
  coordinatesPVDG[0][2]  =  9000000; valuesPVDG_Bg[2]  = 0.01328;  valuesPVDG_visc[2] = 0.00001526;
  coordinatesPVDG[0][3]  = 12000000; valuesPVDG_Bg[3]  = 0.00977;  valuesPVDG_visc[3] = 0.0000166;
  coordinatesPVDG[0][4]  = 15000000; valuesPVDG_Bg[4]  = 0.00773;  valuesPVDG_visc[4] = 0.00001818;
  coordinatesPVDG[0][5]  = 18000000; valuesPVDG_Bg[5]  = 0.006426; valuesPVDG_visc[5] = 0.00001994;
  coordinatesPVDG[0][6]  = 21000000; valuesPVDG_Bg[6]  = 0.005541; valuesPVDG_visc[6] = 0.00002181;
  coordinatesPVDG[0][7]  = 24000000; valuesPVDG_Bg[7]  = 0.004919; valuesPVDG_visc[7] = 0.0000237;
  coordinatesPVDG[0][8]  = 27000000; valuesPVDG_Bg[8]  = 0.004471; valuesPVDG_visc[8] = 0.00002559;
  coordinatesPVDG[0][9]  = 29500000; valuesPVDG_Bg[9]  = 0.004194; valuesPVDG_visc[9] = 0.00002714;
  coordinatesPVDG[0][10] = 31000000; valuesPVDG_Bg[10] = 0.004031; valuesPVDG_visc[10] = 0.00002806;
  coordinatesPVDG[0][11] = 33000000; valuesPVDG_Bg[11] = 0.00391;  valuesPVDG_visc[11] = 0.00002832;
  coordinatesPVDG[0][12] = 53000000; valuesPVDG_Bg[12] = 0.003868; valuesPVDG_visc[12] = 0.00002935;

  TableFunction & tablePVDO_Bo = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "PVDO_Bo" ) );
  tablePVDO_Bo.setTableCoordinates( coordinatesPVDO, { units::Pressure } );
  tablePVDO_Bo.setTableValues( valuesPVDO_Bo, units::Dimensionless );
  tablePVDO_Bo.reInitializeFunction();
  tablePVDO_Bo.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & tablePVDO_visc = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "PVDO_visc" ) );
  tablePVDO_visc.setTableCoordinates( coordinatesPVDO, { units::Pressure } );
  tablePVDO_visc.setTableValues( valuesPVDO_visc, units::Viscosity );
  tablePVDO_visc.reInitializeFunction();
  tablePVDO_visc.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & tablePVDG_Bg = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "PVDG_Bg" ) );
  tablePVDG_Bg.setTableCoordinates( coordinatesPVDG, { units::Pressure } );
  tablePVDG_Bg.setTableValues( valuesPVDG_Bg, units::Dimensionless );
  tablePVDG_Bg.reInitializeFunction();
  tablePVDG_Bg.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & tablePVDG_visc = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "PVDG_visc" ) );
  tablePVDG_visc.setTableCoordinates( coordinatesPVDG, { units::Pressure } );
  tablePVDG_visc.setTableValues( valuesPVDG_visc, units::Viscosity );
  tablePVDG_visc.reInitializeFunction();
  tablePVDG_visc.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then, define the Dead-Oil constitutive model

  DeadOilFluid & fluid = parent->registerGroup< DeadOilFluid >( name );

  string_array & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 3 );
  compNames[0] = "gas"; compNames[1] = "water"; compNames[2] = "oil";

  array1d< real64 > & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 3 );
  molarWgt[0] = 18e-3; molarWgt[1] = 16e-3; molarWgt[2] = 114e-3;

  string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "gas"; phaseNames[1] = "water"; phaseNames[2] = "oil";

  array1d< real64 > & surfaceDens = fluid.getReference< array1d< real64 > >( DeadOilFluid::viewKeyStruct::surfacePhaseMassDensitiesString() );
  surfaceDens.resize( 3 );
  surfaceDens[0] = 0.9907; surfaceDens[1] = 1022.0; surfaceDens[2] = 800.0;

  string_array & FVFTableNames = fluid.getReference< string_array >( DeadOilFluid::viewKeyStruct::formationVolumeFactorTableNamesString() );
  FVFTableNames.resize( 2 );
  FVFTableNames[0] = "PVDG_Bg"; FVFTableNames[1] = "PVDO_Bo";

  string_array & viscosityTableNames = fluid.getReference< string_array >( DeadOilFluid::viewKeyStruct::viscosityTableNamesString() );
  viscosityTableNames.resize( 2 );
  viscosityTableNames[0] = "PVDG_visc"; viscosityTableNames[1] = "PVDO_visc";

  real64 & waterRefPressure = fluid.getReference< real64 >( DeadOilFluid::viewKeyStruct::waterRefPressureString() );
  waterRefPressure = 30600000.1;
  real64 & waterFormationVolumeFactor = fluid.getReference< real64 >( DeadOilFluid::viewKeyStruct::waterFormationVolumeFactorString() );
  waterFormationVolumeFactor = 1.03;
  real64 & waterCompressibility = fluid.getReference< real64 >( DeadOilFluid::viewKeyStruct::waterCompressibilityString() );
  waterCompressibility = 0.00000000041;
  real64 & waterViscosity = fluid.getReference< real64 >( DeadOilFluid::viewKeyStruct::waterViscosityString() );
  waterViscosity = 0.0003;

  fluid.postProcessInputRecursive();
  return fluid;
}

void writeTableToFile( string const & filename, char const * str )
{
  std::ofstream os( filename );
  ASSERT_TRUE( os.is_open() );
  os << str;
  os.close();
}

void removeFile( string const & filename )
{
  int const ret = std::remove( filename.c_str() );
  ASSERT_TRUE( ret == 0 );
}

class LiveOilFluidTest : public CompositionalFluidTestBase
{
public:
  LiveOilFluidTest()
  {
    writeTableToFile( "pvto.txt", pvtoTableContent );
    writeTableToFile( "pvdg.txt", pvdgTableContent );
    writeTableToFile( "pvtw.txt", pvtwTableContent );

    parent.resize( 1 );
    fluid = &makeLiveOilFluid( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~LiveOilFluidTest()
  {
    removeFile( "pvto.txt" );
    removeFile( "pvdg.txt" );
    removeFile( "pvtw.txt" );
  }
};

TEST_F( LiveOilFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  real64 const P[3] = { 5.4e6, 1.24e7, 3.21e7 };
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.79999; comp[1] = 0.2; comp[2] = 0.00001;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-12;

  for( integer i = 0; i < 3; ++i )
  {
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, false, relTol );
  }
}

TEST_F( LiveOilFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  real64 const P[3] = { 5.4e6, 1.24e7, 3.21e7 };
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.79999; comp[1] = 0.2; comp[2] = 0.00001;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-12;

  for( integer i = 0; i < 3; ++i )
  {
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, false, relTol );
  }
}

class DeadOilFluidTest : public CompositionalFluidTestBase
{
public:

  DeadOilFluidTest()
  {
    writeTableToFile( "pvdo.txt", pvdoTableContent );
    writeTableToFile( "pvdg.txt", pvdgTableContent );
    writeTableToFile( "pvdw.txt", pvdwTableContent );

    parent.resize( 1 );
    fluid = &makeDeadOilFluid( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~DeadOilFluidTest()
  {
    removeFile( "pvdo.txt" );
    removeFile( "pvdg.txt" );
    removeFile( "pvdw.txt" );
  }
};

TEST_F( DeadOilFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  real64 const P[3] = { 1.24e7, 3.21e7, 5.01e7 };
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.1; comp[1] = 0.3; comp[2] = 0.6;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;

  for( integer i = 0; i < 3; ++i )
  {
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, false, relTol );
  }
}

TEST_F( DeadOilFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  real64 const P[3] = { 5.4e6, 1.24e7, 3.21e7 };
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.1; comp[1] = 0.3; comp[2] = 0.6;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;
  real64 const absTol = 1e-14;

  for( integer i = 0; i < 3; ++i )
  {
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, false, relTol, absTol );
  }
}

class DeadOilFluidFromTableTest : public CompositionalFluidTestBase
{
public:

  DeadOilFluidFromTableTest()
  {
    parent.resize( 1 );
    fluid = &makeDeadOilFluidFromTable( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }
};

TEST_F( DeadOilFluidFromTableTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  real64 const P[3] = { 5.4e6, 1.24e7, 3.21e7 };
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.1; comp[1] = 0.3; comp[2] = 0.6;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;

  for( integer i = 0; i < 3; ++i )
  {
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, false, relTol );
  }
}

MultiFluidBase & makeCO2BrinePhillipsFluid( string const & name, Group * parent )
{
  CO2BrinePhillipsFluid & fluid = parent->registerGroup< CO2BrinePhillipsFluid >( name );

  auto & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 2 );
  compNames[0] = "co2"; compNames[1] = "water";

  auto & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 2 );
  molarWgt[0] = 44e-3; molarWgt[1] = 18e-3;

  auto & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "gas"; phaseNames[1] = "liquid";

  auto & phasePVTParaFileNames = fluid.getReference< path_array >( CO2BrinePhillipsFluid::viewKeyStruct::phasePVTParaFilesString() );
  phasePVTParaFileNames.resize( 2 );
  phasePVTParaFileNames[0] = "pvtgas.txt"; phasePVTParaFileNames[1] = "pvtliquid.txt";

  auto & flashModelParaFileName = fluid.getReference< Path >( CO2BrinePhillipsFluid::viewKeyStruct::flashModelParaFileString() );
  flashModelParaFileName = "co2flash.txt";

  fluid.postProcessInputRecursive();
  return fluid;
}

class CO2BrinePhillipsFluidTest : public CompositionalFluidTestBase
{
protected:

  CO2BrinePhillipsFluidTest()
  {
    writeTableToFile( "pvtliquid.txt", pvtLiquidPhillipsTableContent );
    writeTableToFile( "pvtgas.txt", pvtGasTableContent );
    writeTableToFile( "co2flash.txt", co2FlashTableContent );

    parent.resize( 1 );
    fluid = &makeCO2BrinePhillipsFluid( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~CO2BrinePhillipsFluidTest()
  {
    removeFile( "pvtliquid.txt" );
    removeFile( "pvtgas.txt" );
    removeFile( "co2flash.txt" );
  }

};


TEST_F( CO2BrinePhillipsFluidTest, checkAgainstPreviousImplementationMolar )
{
  fluid->setMassFlag( false );

  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const T[3] = { 367.65, 368.15, 368.75 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const relTol = 1e-10;

  fluid->allocateConstitutiveData( fluid->getParent(), 1 );

  CO2BrinePhillipsFluid::KernelWrapper wrapper =
    dynamicCast< CO2BrinePhillipsFluid * >( fluid )->createKernelWrapper();

  real64 const savedTotalDens[] =
  { 5881.9010529428224, 5869.6094131788523, 5855.0332090690354, 9181.3523596865525, 9157.6071613646127, 9129.5728206336826, 15757.685798517123, 15698.877814368472,
    15630.149353340639 };
  real64 const savedGasPhaseFrac[] =
  { 0.29413690046142371148, 0.29415754810481165027, 0.29418169867697463449, 0.29194010802017489326, 0.29196434961986583723, 0.29199266189550621142, 0.2890641335638892695, 0.28908718137828937067,
    0.28911404840933618843 };
  real64 const savedWaterDens[] =
  { 53296.719183517576, 53274.578175308554, 53247.856162690216, 53248.577831698305, 53226.801031868054, 53200.505577694363, 53232.959859840405, 53211.345942175059,
    53185.244751356993 };
  real64 const savedGasDens[] =
  { 1876.2436091302606656, 1872.184636376355229, 1867.3711104617746059, 3053.1548401973859654, 3044.5748249030266379, 3034.4507978134674886, 5769.0622621289458039, 5742.8476745352018042,
    5712.2837704249559465 };
  real64 const savedWaterMassDens[] =
  { 970.85108546544745423, 970.4075834766143771, 969.87385780866463847, 974.23383396044232541, 973.78856424100911227, 973.25280170872576946, 979.48333010951580491, 979.04147229150635212,
    978.50977403260912979 };
  real64 const savedGasMassDens[] =
  { 82.554718801731468147, 82.376124000559627802, 82.164328860318079251, 134.33881296868497657, 133.96129229573315911, 133.51583510379256836, 253.83873953367358922, 252.68529767954885301,
    251.34048589869803436 };
  real64 const savedWaterVisc[] =
  { 0.0003032144206279845924, 0.00030157070452334377216, 0.00029959815370251820189, 0.0003032144206279845924, 0.00030157070452334377216, 0.00029959815370251820189, 0.0003032144206279845924,
    0.00030157070452334377216, 0.00029959815370251820189 };

  real64 const savedGasVisc[] =
  { 1.9042384704865343673e-05, 1.9062615947696152414e-05, 1.9086923154230274463e-05, 2.0061713844617985449e-05, 2.0075955757102255573e-05, 2.0093249989250199265e-05, 2.3889596884008691474e-05,
    2.3865756080512667728e-05, 2.3839170076324036522e-05  };
  real64 const savedWaterPhaseGasComp[] =
  { 0.0083062842389820552153, 0.008277274736736653718, 0.0082433415400525456018, 0.011383065290266058955, 0.011349217198060387521, 0.011309682362800700661, 0.015382352969377973903,
    0.015350431636424789056, 0.015313218057419366105  };
  real64 const savedWaterPhaseWaterComp[] =
  { 0.99169371576101794652, 0.9917227252632633272, 0.99175665845994742664, 0.98861693470973388553, 0.98865078280193963156, 0.98869031763719927852, 0.98461764703062204518, 0.98464956836357520054,
    0.98468678194258063563 };

  integer counter = 0;
  for( integer i = 0; i < 3; ++i )
  {
    for( integer j = 0; j < 3; ++j )
    {
      testValuesAgainstPreviousImplementation( wrapper,
                                               P[i], T[j], comp,
                                               savedTotalDens[counter], savedGasPhaseFrac[counter],
                                               savedWaterDens[counter], savedGasDens[counter],
                                               savedWaterMassDens[counter], savedGasMassDens[counter],
                                               savedWaterVisc[counter], savedGasVisc[counter],
                                               savedWaterPhaseGasComp[counter], savedWaterPhaseWaterComp[counter],
                                               relTol );
      counter++;
    }
  }
}

TEST_F( CO2BrinePhillipsFluidTest, checkAgainstPreviousImplementationMass )
{
  fluid->setMassFlag( true );

  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const T[3] = { 367.65, 368.15, 368.75 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const relTol = 1e-10;

  fluid->allocateConstitutiveData( fluid->getParent(), 1 );

  CO2BrinePhillipsFluid::KernelWrapper wrapper =
    dynamicCast< CO2BrinePhillipsFluid * >( fluid )->createKernelWrapper();

  real64 const savedTotalDens[] =
  { 238.31504112633914, 237.83897306400553, 237.27445306546298, 353.95258514794097, 353.12773295711992, 352.1532278769692, 549.90502586392017, 548.2725957521294,
    546.35992000222234 };
  real64 const savedGasPhaseFrac[] =
  { 0.28566797890570228, 0.28571845092287312, 0.28577748565482669, 0.28029804182709406, 0.28035729907078311, 0.2804265068556816, 0.27326788204506247, 0.27332422114692956,
    0.27338989611171055 };
  real64 const savedWaterDens[] =
  { 970.85108546544745423, 970.4075834766143771, 969.87385780866463847, 974.23383396044232541, 973.78856424100911227, 973.25280170872576946, 979.48333010951580491, 979.04147229150635212,
    978.50977403260912979 };
  real64 const savedGasDens[] =
  { 82.554718801731468147, 82.376124000559627802, 82.164328860318079251, 134.33881296868497657, 133.96129229573315911, 133.51583510379256836, 253.83873953367358922, 252.68529767954885301,
    251.34048589869803436 };
  real64 const savedWaterMassDens[] =
  { 970.85108546544745423, 970.4075834766143771, 969.87385780866463847, 974.23383396044232541, 973.78856424100911227, 973.25280170872576946, 979.48333010951580491, 979.04147229150635212,
    978.50977403260912979 };
  real64 const savedGasMassDens[] =
  { 82.554718801731468147, 82.376124000559627802, 82.164328860318079251, 134.33881296868497657, 133.96129229573315911, 133.51583510379256836, 253.83873953367358922, 252.68529767954885301,
    251.34048589869803436 };
  real64 const savedWaterVisc[] =
  { 0.0003032144206279845924, 0.00030157070452334377216, 0.00029959815370251820189, 0.0003032144206279845924, 0.00030157070452334377216, 0.00029959815370251820189, 0.0003032144206279845924,
    0.00030157070452334377216, 0.00029959815370251820189 };
  real64 const savedGasVisc[] =
  { 1.9042384704865343673e-05, 1.9062615947696152414e-05, 1.9086923154230274463e-05, 2.0061713844617985449e-05, 2.0075955757102255573e-05, 2.0093249989250199265e-05, 2.3889596884008691474e-05,
    2.3865756080512667728e-05, 2.3839170076324036522e-05  };
  real64 const savedWaterPhaseGasComp[] =
  { 0.020063528822832473, 0.019994285300494085, 0.019913282008777046, 0.027375162661670061, 0.027295074213708581, 0.02720152052681666, 0.036784005129927563,
    0.036709327088310859, 0.036622259649145526 };
  real64 const savedWaterPhaseWaterComp[] =
  { 0.97993647117716742, 0.98000571469950604, 0.98008671799122282, 0.97262483733832983, 0.97270492578629131, 0.97279847947318321, 0.96321599487007259, 0.96329067291168902,
    0.96337774035085455 };

  integer counter = 0;
  for( integer i = 0; i < 3; ++i )
  {
    for( integer j = 0; j < 3; ++j )
    {
      testValuesAgainstPreviousImplementation( wrapper,
                                               P[i], T[j], comp,
                                               savedTotalDens[counter], savedGasPhaseFrac[counter],
                                               savedWaterDens[counter], savedGasDens[counter],
                                               savedWaterMassDens[counter], savedGasMassDens[counter],
                                               savedWaterVisc[counter], savedGasVisc[counter],
                                               savedWaterPhaseGasComp[counter], savedWaterPhaseWaterComp[counter],
                                               relTol );
      counter++;
    }
  }
}

TEST_F( CO2BrinePhillipsFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  // TODO test over a range of values
  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const T[3] = { 367.65, 368.15, 368.75 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;

  for( integer i = 0; i < 3; ++i )
  {
    for( integer j = 0; j < 3; ++j )
    {
      testNumericalDerivatives( *fluid, parent, P[i], T[j], comp, eps, false, relTol );
    }
  }
}

TEST_F( CO2BrinePhillipsFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  // TODO test over a range of values
  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const T[3] = { 367.65, 368.15, 368.75 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-8;

  for( integer i = 0; i < 3; ++i )
  {
    for( integer j = 0; j < 3; ++j )
    {
      testNumericalDerivatives( *fluid, parent, P[i], T[j], comp, eps, false, relTol );
    }
  }
}

MultiFluidBase & makeCO2BrineEzrokhiFluid( string const & name, Group * parent )
{
  CO2BrineEzrokhiFluid & fluid = parent->registerGroup< CO2BrineEzrokhiFluid >( name );

  auto & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 2 );
  compNames[0] = "co2"; compNames[1] = "water";

  auto & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 2 );
  molarWgt[0] = 44e-3; molarWgt[1] = 18e-3;

  auto & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "gas"; phaseNames[1] = "liquid";

  auto & phasePVTParaFileNames = fluid.getReference< path_array >( CO2BrineEzrokhiFluid::viewKeyStruct::phasePVTParaFilesString() );
  phasePVTParaFileNames.resize( 2 );
  phasePVTParaFileNames[0] = "pvtgas.txt"; phasePVTParaFileNames[1] = "pvtliquid.txt";

  auto & flashModelParaFileName = fluid.getReference< Path >( CO2BrineEzrokhiFluid::viewKeyStruct::flashModelParaFileString() );
  flashModelParaFileName = "co2flash.txt";

  fluid.postProcessInputRecursive();
  return fluid;
}

class CO2BrineEzrokhiFluidTest : public CompositionalFluidTestBase
{
protected:

  CO2BrineEzrokhiFluidTest()
  {
    writeTableToFile( "pvtliquid.txt", pvtLiquidEzrokhiTableContent );
    writeTableToFile( "pvtgas.txt", pvtGasTableContent );
    writeTableToFile( "co2flash.txt", co2FlashTableContent );

    parent.resize( 1 );
    fluid = &makeCO2BrineEzrokhiFluid( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~CO2BrineEzrokhiFluidTest()
  {
    removeFile( "pvtliquid.txt" );
    removeFile( "pvtgas.txt" );
    removeFile( "co2flash.txt" );
  }

};

TEST_F( CO2BrineEzrokhiFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  // TODO test over a range of values
  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const T[3] = { 367.65, 368.15, 368.75 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;

  for( integer i = 0; i < 3; ++i )
  {
    for( integer j = 0; j < 3; ++j )
    {
      testNumericalDerivatives( *fluid, parent, P[i], T[j], comp, eps, false, relTol );
    }
  }
}

TEST_F( CO2BrineEzrokhiFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  // TODO test over a range of values
  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const T[3] = { 367.65, 368.15, 368.75 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-8;

  for( integer i = 0; i < 3; ++i )
  {
    for( integer j = 0; j < 3; ++j )
    {
      testNumericalDerivatives( *fluid, parent, P[i], T[j], comp, eps, false, relTol );
    }
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
