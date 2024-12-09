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

#include "constitutiveTestHelpers.hpp"

#include "constitutive/fluid/multifluid/compositional/models/KValueFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/KValueFlashParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/models/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"
#include "constitutive/unitTests/TestFluid.hpp"
#include "constitutive/unitTests/TestFluidUtilities.hpp"

#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geos::dataRepository;
using namespace geos::testing;
using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{

template< integer numPhases, integer numComps >
using FlashData = std::tuple<
  real64 const,                 // pressure
  real64 const,                 // temperature
  Feed< numComps > const,       // total composition
  Feed< numPhases > const,      // expected phase fractions
  Feed< 2*numPhases >           // expected phase composition (2 selected components)
  >;

template< int numComps >
struct TestData;

template<>
struct TestData< 9 >
{
  static std::unique_ptr< TestFluid< 9 > > createFluid()
  {
    auto fluid = TestFluid< 9 >::create( {Fluid::H2O, Fluid::CO2, Fluid::N2, Fluid::C1, Fluid::C2, Fluid::C3, Fluid::C4, Fluid::C5, Fluid::C10} );
    const std::array< real64 const, 36 > bics = {
      0.01, 0, 0.003732, 0, 0.01, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0.01, 0, 0.028, 0.01, 0.01, 0, 0, 0.01, 0, 0.04532, 0.01, 0.01, 0, 0, 0
    };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template< integer numPhases, integer numComps >
class KValueFlashTestFixture : public ConstitutiveTestBase< MultiFluidBase >, public ::testing::WithParamInterface< FlashData< numPhases, numComps > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numDofs = numComps + 2;
  // Selected components for test
  static constexpr int comp0 = 0;
  static constexpr int comp1 = numComps-1;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
  using FlashModelType = KValueFlashModel< numPhases >;
  using FlashModelParamType = KValueFlashParameters< numPhases >;
  using PhasePropSlice = typename FlashModelType::KernelWrapper::PhaseProp::SliceType;
  using PhaseCompSlice = typename FlashModelType::KernelWrapper::PhaseComp::SliceType;
public:
  KValueFlashTestFixture();
  ~KValueFlashTestFixture() override;

  void testFlash( typename KValueFlashTestFixture::ParamType const & data );

protected:
  std::unique_ptr< TestFluid< numComps > > m_fluid{};
  std::unique_ptr< ModelParameters > m_parameters{};
  std::unique_ptr< FlashModelType > m_flash{};

private:
  void generateTables( arraySlice1d< string > const & names, string const fluidName );

  static void writeToFile( string const & fileName, string const & content );

  static void removeFile( string const & fileName );

  static CompositionalKValueConstantViscosity * makeFluid( string const & name,
                                                           Group * parent,
                                                           TestFluid< numComps > const * testFluid,
                                                           FlashModelParamType const * parameters );

private:
  string_array m_fileNames;
};

template< integer numPhases, integer numComps >
KValueFlashTestFixture< numPhases, numComps >::KValueFlashTestFixture():
  m_fluid( TestData< numComps >::createFluid() )
{
  string const fluidName = GEOS_FMT( "fluid_{}_{}", numPhases, numComps );

  m_parameters = FlashModelType::createParameters( std::move( m_parameters ));
  FlashModelParamType * parameters = const_cast< FlashModelParamType * >(m_parameters->get< FlashModelParamType >());
  parameters->m_kValueTables.resize( (numPhases-1)*numComps );
  generateTables( parameters->m_kValueTables.toSlice(), fluidName );

  ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

  auto & parent = this->m_parent;
  parent.resize( 1 );

  m_model = makeFluid( fluidName, &parent, m_fluid.get(), parameters );

  parent.initialize();
  parent.initializePostInitialConditions();

  m_parameters->postInputInitialization( m_model, componentProperties );

  string const flashName = GEOS_FMT( "{}_flash_copy", fluidName );
  m_flash = std::make_unique< FlashModelType >( flashName, componentProperties, *m_parameters );
}

template< integer numPhases, integer numComps >
KValueFlashTestFixture< numPhases, numComps >::~KValueFlashTestFixture()
{
  for( string const & fileName : m_fileNames )
  {
    removeFile( fileName );
  }
}

template< integer numComps >
struct MakeFluid;

template<>
struct MakeFluid< 9 >
{
  static void populate( CompositionalKValueConstantViscosity & fluid, TestFluid< 9 > const * testFluid )
  {
    using FluidModel = CompositionalKValueConstantViscosity;

    string_array & componentNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    TestFluid< 9 >::createArray( componentNames, testFluid->componentNames );

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    TestFluid< 9 >::createArray( molarWeight, testFluid->molecularWeight );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentCriticalPressureString() );
    TestFluid< 9 >::createArray( criticalPressure, testFluid->criticalPressure );

    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentCriticalTemperatureString() );
    TestFluid< 9 >::createArray( criticalTemperature, testFluid->criticalTemperature );

    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentAcentricFactorString() );
    TestFluid< 9 >::createArray( acentricFactor, testFluid->acentricFactor );
  }
};

template< integer numPhases, integer numComps >
CompositionalKValueConstantViscosity *
KValueFlashTestFixture< numPhases, numComps >::makeFluid( string const & name,
                                                          Group * parent,
                                                          TestFluid< numComps > const * testFluid,
                                                          FlashModelParamType const * parameters )
{
  CompositionalKValueConstantViscosity & compositionalFluid = parent->registerGroup< CompositionalKValueConstantViscosity >( name );

  Group & fluid = compositionalFluid;

  auto & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.emplace_back( "gas" );
  phaseNames.emplace_back( "liquid" );

  string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
  string_array & equationsOfState = fluid.template getReference< string_array >( EquationOfState::viewKeyStruct::equationsOfStateString() );
  equationsOfState.emplace_back( eosName );
  equationsOfState.emplace_back( eosName );

  MakeFluid< numComps >::populate( compositionalFluid, testFluid );

  string_array & kValueTables = fluid.template getReference< string_array >( FlashModelParamType::viewKeyStruct::kValueTablesString() );
  for( auto const & tableName : parameters->m_kValueTables )
  {
    kValueTables.emplace_back( tableName );
  }

  compositionalFluid.postInputInitializationRecursive();

  return &compositionalFluid;
}

// Crookston correlations pressure (bar), temperature (K)
real64 getKValue( integer const phaseIndex, integer const compIndex, real64 const pressure, real64 const temperature )
{
  GEOS_UNUSED_VAR( phaseIndex );

  static std::array< real64, 5 > constexpr crookstonCoefficients[9] = {
    {3.0620e+00, 8.9414e+02, 1.1912e-02, 5.3659e+02, 1.1951e+02},
    {-1.8141e+00, 6.2655e+02, 6.7489e-03, 5.0732e+00, 2.5249e+02},
    {-3.5742e-01, 4.8660e+02, 5.7887e-03, 9.1910e+01, 1.8027e+02},
    {3.6577e+00, 7.1889e+02, 1.4154e-02, 5.7323e+02, 1.1088e+02},
    {8.2845e+00, 1.0144e+03, 5.2968e-02, 9.5800e+02, 9.3548e+01},
    {1.5732e+01, 1.4968e+03, 1.7407e-01, 1.3707e+03, 8.1261e+01},
    {2.5068e+01, 2.2625e+03, 5.3888e-01, 1.7981e+03, 7.1948e+01},
    {3.6243e+01, 4.1985e+03, 1.9055e+00, 2.3575e+03, 5.8025e+01},
    {-2.3346e-01, 7.9356e-01, 6.8406e-03, 7.1217e+02, 1.3827e+02}
  };
  auto const [A, B, C, D, E] = std::apply( []( auto &&... args ) { return std::make_tuple( args ... ); }, crookstonCoefficients[compIndex] );
  return (A + B/pressure + C*pressure)*LvArray::math::exp( -D/(temperature-E));
}

template< integer numPhases, integer numComps >
void KValueFlashTestFixture< numPhases, numComps >::generateTables( arraySlice1d< string > const & names, string const fluidName )
{
  FunctionManager & functionManager = FunctionManager::getInstance();
  std::ostringstream content;
  for( integer ip = 0; ip < numPhases-1; ++ip )
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      string const tableName = GEOS_FMT( "{}_KVALUE_{}_{}", fluidName, ip+1, ic+1 );

      // Generate tables whose end-points are not exactly the same
      real64 const th = (ip + 0.5)*numComps + ic + 0.5;
      real64 const c = LvArray::math::cos( th );
      real64 const s = LvArray::math::sin( th );
      // Pressure in bar
      real64 const minPressure = 1.0 + c*c;
      real64 const maxPressure = 600.0 + 50.0*c*c;
      integer const NP = static_cast< integer >(6.0 + 8.0*s*s);
      real64 const dp = (maxPressure - minPressure)/NP;
      content.str( "" );
      for( integer i = 0; i <= NP; ++i )
      {
        real64 const pressure = 1.0e5 * (minPressure + i*dp);
        content << pressure << "\n";
      }
      string const pressureFileName = GEOS_FMT( "{}_PRESSURE.txt", tableName );
      writeToFile( pressureFileName, content.str());
      m_fileNames.emplace_back( pressureFileName );

      // Temp in degC
      real64 const minTemp = 10.0 + 5.0*s*s;
      real64 const maxTemp = 200.0 + 50.0*s*s;
      integer const NT = static_cast< integer >(6.0 + 8.0*c*c);
      real64 const dt = (maxTemp - minTemp)/NP;
      content.str( "" );
      for( integer j = 0; j <= NT; ++j )
      {
        real64 const temperature = minTemp + j*dt + 273.15;
        content << temperature << "\n";
      }
      string const temperatureFileName = GEOS_FMT( "{}_TEMP.txt", tableName );
      writeToFile( temperatureFileName, content.str());
      m_fileNames.emplace_back( temperatureFileName );

      content.str( "" );
      for( integer i = 0; i <= NP; ++i )
      {
        real64 const pressure = minPressure + i*dp;
        for( integer j = 0; j <= NT; ++j )
        {
          real64 const temperature = minTemp + j*dt + 273.15;
          real64 const kValue = getKValue( ip, ic, pressure, temperature );
          content << kValue << "\n";
        }
      }
      string const kValueFileName = GEOS_FMT( "{}_KVALUE.txt", tableName );
      writeToFile( kValueFileName, content.str());
      m_fileNames.emplace_back( kValueFileName );

      TableFunction * tableFunction = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", tableName ) );

      path_array & coordinateFiles = tableFunction->getWrapper< path_array >( TableFunction::viewKeyStruct::coordinateFilesString())
                                       .reference();
      coordinateFiles.emplace_back( Path() );
      coordinateFiles.emplace_back( Path() );
      dynamicCast< string & >( coordinateFiles[0] ) = pressureFileName;
      dynamicCast< string & >( coordinateFiles[1] ) = temperatureFileName;
      Path & voxelFile = tableFunction->getWrapper< Path >( TableFunction::viewKeyStruct::voxelFileString())
                           .reference();
      dynamicCast< string & >( voxelFile ) = kValueFileName;
      tableFunction->setInterpolationMethod( TableFunction::InterpolationType::Linear );
      tableFunction->initializeFunction();

      integer const tableIndex = numComps*ip + ic;

      names[tableIndex] = tableName;
    }
  }
}

template< integer numPhases, integer numComps >
void KValueFlashTestFixture< numPhases, numComps >::writeToFile( string const & fileName, string const & content )
{
  std::ofstream os( fileName );
  ASSERT_TRUE( os.is_open() );
  os << content;
  os.close();
}

template< integer numPhases, integer numComps >
void KValueFlashTestFixture< numPhases, numComps >::removeFile( string const & fileName )
{
  int const ret = std::remove( fileName.c_str() );
  ASSERT_TRUE( ret == 0 );
}

/* --- Start tests --- */

template< integer numPhases, integer numComps >
void KValueFlashTestFixture< numPhases, numComps >::testFlash( typename KValueFlashTestFixture::ParamType const & data )
{
  ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

  auto flashKernelWrapper = this->m_flash->createKernelWrapper();
  GEOS_UNUSED_VAR( flashKernelWrapper );

  real64 const pressure = std::get< 0 >( data );
  real64 const temperature = std::get< 1 >( data );
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, std::get< 2 >( data ));

  stackArray2d< real64, (numPhases-1)*numComps > kValues( numPhases-1, numComps );
  LvArray::forValuesInSlice( kValues.toSlice(), []( real64 & v ){ v = 0.0; } );

  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseFractionData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases *numDofs, multifluid::LAYOUT_PHASE_DC > dPhaseFractionData( 1, 1, numPhases, numDofs );
  StackArray< real64, 4, numPhases *numComps, multifluid::LAYOUT_PHASE_COMP > phaseComponentFractionData( 1, 1, numPhases, numComps );
  StackArray< real64, 5, numPhases *numComps *numDofs, multifluid::LAYOUT_PHASE_COMP_DC > dPhaseComponentFractionData( 1, 1, numPhases, numComps, numDofs );

  auto phaseFraction = phaseFractionData[0][0];
  auto dPhaseFraction = dPhaseFractionData[0][0];
  auto phaseComponentFraction = phaseComponentFractionData[0][0];
  auto dPhaseComponentFraction = dPhaseComponentFractionData[0][0];

  flashKernelWrapper.compute( componentProperties.createKernelWrapper(),
                              pressure,
                              temperature,
                              composition.toSliceConst(),
                              kValues.toSlice(),
                              PhasePropSlice( phaseFraction, dPhaseFraction ),
                              PhaseCompSlice( phaseComponentFraction, dPhaseComponentFraction ) );

}

using KValueFlashTest_2_9 = KValueFlashTestFixture< 2, 9 >;
TEST_P( KValueFlashTest_2_9, testFlash )
{
  testFlash( GetParam() );
}

/* UNCRUSTIFY-OFF */

// Test data
INSTANTIATE_TEST_SUITE_P(
  KValueFlash, KValueFlashTest_2_9,
  ::testing::Values( 
    FlashData<2, 9>( 1.0e+06, 298.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.000197, 0.999803}, {0.000361, 0.839175, 0.010852, 0.000016} ),
    FlashData<2, 9>( 1.0e+07, 298.15, {0.000000, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.999999}, {0.000197, 0.999803}, {0.000361, 0.839175, 0.010852, 0.000016} ),
    FlashData<2, 9>( 1.0e+06, 298.15, {0.100000, 0.000000, 0.899999, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001}, {0.000197, 0.999803}, {0.000361, 0.839175, 0.010852, 0.000016} ),
    FlashData<2, 9>( 1.0e+07, 298.15, {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000}, {0.000197, 0.999803}, {0.000361, 0.839175, 0.010852, 0.000016} ),
    FlashData<2, 9>( 1.0e+06, 298.15, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, {0.000197, 0.999803}, {0.000361, 0.839175, 0.010852, 0.000016} )
  )
);

/* UNCRUSTIFY-ON */

} // namespace geos

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
