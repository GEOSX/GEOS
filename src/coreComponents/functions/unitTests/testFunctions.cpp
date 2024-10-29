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

#include "codingUtilities/UnitTestUtilities.hpp"
#include "gtest/gtest.h"
#include "mainInterface/initialization.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/FunctionBase.hpp"
#include "functions/TableFunction.hpp"
#include "functions/MultivariableTableFunction.hpp"
#include "functions/MultivariableTableFunctionKernels.hpp"
//#include "mainInterface/GeosxState.hpp"

#ifdef GEOS_USE_MATHPRESSO
  #include "functions/SymbolicFunction.hpp"
#endif

#include <random>

using namespace geos;

char const * multivariableTableFileContent =
  "4 34\n"
  "3 1 1000\n"
  "2 1e-13 0.99999999999989997\n"
  "2 1e-13 0.99999999999989997\n"
  "2 1e-13 0.99999999999989997\n"
  "2.5165096633488189e-12 2.5165096633488189e-12 2.5165096633488189e-12 19.982814779283817 0.17715834234447858 3.4799790238502815 3.4638346743277832 9.9999999999999998e-13 46.589817605162132 2.28794740689822e-11 0.22773331410153216 9.9999999999999998e-13 1.2444262762971973e-13 1.7553963362399752e-13 0.00044473177341675438 0.0058778098301284758 0 0 0 0 0 0 0 -4.5296539991171088 0 2.2648269995585544 0 0 0 1419.9027145548464 600 0 0 2.999822612537173e-13 \n"
  "7.7713016930315248e-12 7.7713016930315248e-12 77.713016930307475 -1.9966842958308671e-12 7.7713016930307477e-12 7.7713016930307477e-12 77.713016930299716 9.9999999999999998e-13 0 0 0 0 1.0000000000000999 0 3.3572023313892832e-15 0 0 0 0 0 0 0 0 0 0 -0 0 0 0 1400.00000000008 0 0 0 1.0000000000000999 \n"
  "2.7975940691009456e-12 27.975940691006656 2.7975940691009456e-12 -1.9966842958308671e-12 2.7975940691006661e-12 27.97594069100386 2.7975940691006661e-12 9.9999999999999998e-13 0 0 0 0 1.0000000000000999 0 1.2085606378514876e-15 0 0 0 0 0 0 0 0 0 0 -0 0 0 0 1400.00000000008 0 0 0 1.0000000000000999 \n"
  "4.1141379411678857e-12 41.141379411674741 41.141379411674741 -19.982814779287814 2.0570689705840455e-12 20.570689705838397 20.570689705838397 9.9999999999999998e-13 0 0 0 0 1.9999999999999001 0 8.8865379529230774e-16 0 0 0 0 0 0 0 0 0 0 -0 0 0 0 1400.00000000004 0 0 0 1.9999999999999001 \n"
  "13.633265167006751 1.3633265167008117e-12 1.3633265167008117e-12 -1.9966842958308671e-12 0 0 0 0 136.33265167005388 1.3633265167006753e-11 1.3633265167006753e-11 9.9999999999999998e-13 0 1.0000000000000999 0 0.005889570552146328 0 0 0 0 0 0 0 -4.9960036108132044e-13 0 2.4980018054066022e-13 0 0 0 0 600 0 0 1.0000000000000999 \n"
  "23.490714960798009 2.3490714960800358e-12 23.490714960798009 -19.982814779287814 0.046287023642824505 1.8411044362036712e-13 1.8234282041118481 9.9999999999999998e-13 96.552512027963189 9.6011389494930567e-24 0.95089595179084896 9.9999999999999998e-13 0.31350334619913905 1.6864966538007611 0.00081380177433002497 0.0058659228652333875 0 0 0 0 0 0 0 -4.9999999999995008 0 2.4999999999997504 0 0 0 1419.8049512378052 600 0 0 1.9999999999999001 \n"
  "18.603704313963757 18.603704313963757 1.8603704313965617e-12 -19.982814779287811 0.07994912869223153 3.1180160190001471 3.0886736196170012e-13 9.9999999999999998e-13 60.252762370057184 5.874644331086449e-11 5.8193604072238805e-14 9.9999999999999998e-13 0.67040756322114736 1.3295924367787526 0.00030738287116134427 0.0058895705521407066 0 0 0 0 0 0 0 33.024999999975314 0 -16.512499999987657 0 0 0 1419.9999999999795 600 0 0 1.9999999999998999 \n"
  "25.166584712355728 25.166584712355728 25.166584712355728 -39.965629558573625 0.17715834234447866 3.4799790238502832 3.463834674327785 9.9999999999999998e-13 46.58981760516216 2.2879474068982213e-11 0.22773331410153225 9.9999999999999998e-13 1.2444998625214401 1.7555001374782597 0.00044473177341675433 0.0058778098301284767 0 0 0 0 0 0 0 9.0593079982355817 0 -4.5296539991177909 0 0 0 1419.9027145548464 600 0 0 2.9999999999996998 \n"
  "2.5892526305777112e-12 2.5892526305777112e-12 2.5892526305777112e-12 19.983812920882041 0.18743599438218472 3.6818662905312016 3.6647853438122544 9.9999999999999998e-13 46.971091776178703 2.3066711387292285e-11 0.22959700095437138 9.9999999999999998e-13 1.2797587881148915e-13 1.7202136655617776e-13 0.00044495391693757605 0.0061714064311433935 0 0 0 0 0 0 0 -4.5296539991171088 0 2.2648269995585544 0 0 0 1420.6119559607666 629.96999999999991 0 0 2.999822612537173e-13 \n"
  "7.775571828640934e-12 7.775571828640934e-12 77.75571828640156 -1.9967840302114436e-12 7.775183458226419e-12 7.775183458226419e-12 77.751834582256407 9.9999999999999998e-13 0 0 0 0 1.0000499500000999 0 3.3588792539538123e-15 0 0 0 0 0 0 0 0 0 0 -0 0 0 0 1400.6993000000803 0 0 0 1.0000000000000999 \n"
  "2.7991312769622549e-12 27.991312769619753 2.7991312769622549e-12 -1.9967840302114436e-12 2.7989914673381822e-12 27.989914673379023 2.7989914673381822e-12 9.9999999999999998e-13 0 0 0 0 1.0000499500000999 0 1.2091643138900945e-15 0 0 0 0 0 0 0 0 0 0 -0 0 0 0 1400.6993000000803 0 0 0 1.0000000000000999 \n"
  "4.1163985569075047e-12 41.163985569070931 41.163985569070931 -19.983812920886038 2.0580964765348525e-12 20.580964765346465 20.580964765346465 9.9999999999999998e-13 0 0 0 0 2.0000998999999 0 8.8909767786305609e-16 0 0 0 0 0 0 0 0 0 0 -0 0 0 0 1400.6993000000402 0 0 0 1.9999999999999001 \n"
  "14.314961758724504 1.4314961758725935e-12 1.4314961758725935e-12 -1.9967840302114436e-12 0 0 0 0 143.14246762097306 1.4314246762098739e-11 1.4314246762098739e-11 9.9999999999999998e-13 0 1.0000499500000999 0 0.0061837546012260356 0 0 0 0 0 0 0 -4.9960036108132044e-13 0 2.4980018054066022e-13 0 0 0 0 629.96999999999991 0 0 1.0000000000000999 \n"
  "24.475680879962209 2.4475680879964658e-12 24.475680879962209 -19.983812920886038 0.050219928346927054 1.9975389556008147e-13 1.9783608136674524 9.9999999999999998e-13 99.82255200313746 9.9263102734970568e-24 0.98310089094030273 9.9999999999999998e-13 0.32648546522983846 1.6736144347700614 0.00081420826831630274 0.0061589257123517939 0 0 0 0 0 0 0 -4.9999999999995008 0 2.4999999999997504 0 0 0 1420.5141438109486 629.96999999999991 0 0 1.9999999999999001 \n"
  "19.215576739407531 19.215576739407531 1.9215576739409455e-12 -19.983812920886034 0.085243537804597286 3.3244979743826186 3.2932124560535381e-13 9.9999999999999998e-13 61.217128048687037 5.9686699847529531e-11 5.9125012231392857e-14 9.9999999999999998e-13 0.69211143544999387 1.3079884645499058 0.0003075364089054894 0.0061837546012201349 0 0 0 0 0 0 0 33.024999999975314 0 -16.512499999987657 0 0 0 1420.7092899999798 629.96999999999991 0 0 1.9999999999998999 \n"
  "25.894057399489316 25.894057399489316 25.894057399489316 -39.967625841770072 0.18743599438218464 3.6818662905312003 3.664785343812254 9.9999999999999998e-13 46.971091776178696 2.3066711387292279e-11 0.2295970009543713 9.9999999999999998e-13 1.2798344636442114 1.7203153863554879 0.00044495391693757605 0.0061714064311433935 0 0 0 0 0 0 0 9.0593079982355817 0 -4.5296539991177909 0 0 0 1420.6119559607666 629.96999999999991 0 0 2.9999999999996998 \n"
  "2.6591944680280507e-12 2.6591944680280507e-12 2.6591944680280507e-12 19.984811062480269 0.19758054956862175 3.8811390924094624 3.8631336775423457 9.9999999999999998e-13 47.288336081028589 2.3222504717668852e-11 0.23114770668014989 9.9999999999999998e-13 1.3136722480775809e-13 1.6864500467385845e-13 0.00044517606045839766 0.0064650030321583112 0 0 0 0 0 0 0 -4.5296539991171088 0 2.2648269995585544 0 0 0 1421.3211973666866 659.94000000000005 0 0 2.999822612537173e-13 \n"
  "7.779842352038681e-12 7.779842352038681e-12 77.798423520379032 -1.9968837645920206e-12 7.7790652234220855e-12 7.7790652234220855e-12 77.79065223421307 9.9999999999999998e-13 0 0 0 0 1.0000999000000999 0 3.3605561765183409e-15 0 0 0 0 0 0 0 0 0 0 -0 0 0 0 1401.39860000008 0 0 0 1.0000000000000999 \n"
  "2.8006686244236482e-12 28.006686244233684 2.8006686244236482e-12 -1.9968837645920206e-12 2.8003888655756975e-12 28.003888655754174 2.8003888655756975e-12 9.9999999999999998e-13 0 0 0 0 1.0000999000000999 0 1.2097679899287013e-15 0 0 0 0 0 0 0 0 0 0 -0 0 0 0 1401.39860000008 0 0 0 1.0000000000000999 \n"
  "4.1186593779428123e-12 41.186593779424001 41.186593779424001 -19.984811062484265 2.0591239824856591e-12 20.591239824854529 20.591239824854529 9.9999999999999998e-13 0 0 0 0 2.0001997999998999 0 8.8954156043380464e-16 0 0 0 0 0 0 0 0 0 0 -0 0 0 0 1401.39860000004 0 0 0 1.9999999999999001 \n"
  "14.99672638050361 1.4996726380505108e-12 1.4996726380505108e-12 -1.9968837645920206e-12 0 0 0 0 149.95228357189228 1.4995228357190726e-11 1.4995228357190726e-11 9.9999999999999998e-13 0 1.0000999000000999 0 0.0064779386503057467 0 0 0 0 0 0 0 -4.9960036108132044e-13 0 2.4980018054066022e-13 0 0 0 0 659.94000000000005 0 0 1.0000000000000999 \n"
  "25.445926657502518 2.5445926657505065e-12 25.445926657502518 -19.984811062484265 0.054247898177082414 2.157754768577228e-13 2.1370384130372595 9.9999999999999998e-13 102.983507562632 1.0240634291613648e-23 1.0142315138747284 9.9999999999999998e-13 0.33925837145800958 1.6609414285418904 0.00081461476230258061 0.0064519285594702028 0 0 0 0 0 0 0 -4.9999999999995008 0 2.4999999999997504 0 0 0 1421.2233363840917 659.94000000000005 0 0 1.9999999999999001 \n"
  "19.80823163954285 19.80823163954285 1.980823163954483e-12 -19.984811062484262 0.090528615667852308 3.5306160110497706 3.4973907984679659e-13 9.9999999999999998e-13 62.091127816030024 6.0538849620689803e-11 5.9969142764489583e-14 9.9999999999999998e-13 0.7131018109295213 1.2870979890703782 0.00030768994664963448 0.0064779386502995641 0 0 0 0 0 0 0 33.024999999975314 0 -16.512499999987657 0 0 0 1421.4185799999796 659.94000000000005 0 0 1.9999999999998999 \n"
  "26.593517132455098 26.593517132455098 26.593517132455098 -39.969622124966527 0.19758054956862175 3.8811390924094633 3.8631336775423457 9.9999999999999998e-13 47.288336081028568 2.3222504717668846e-11 0.23114770668014975 9.9999999999999998e-13 1.313749928999681 1.6865497710000181 0.0004451760604583976 0.0064650030321583112 0 0 0 0 0 0 0 9.0593079982355817 0 -4.5296539991177909 0 0 0 1421.3211973666866 659.94000000000005 0 0 2.9999999999996998";



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

void evaluate1DFunction( FunctionBase & function,
                         arrayView1d< real64 const > const & inputs,
                         arrayView1d< real64 const > const & outputs )
{
  for( localIndex ii=0; ii<inputs.size(); ++ii )
  {
    real64 input = inputs[ii];
    real64 predicted = function.evaluate( &input );
    real64 expected = outputs[ii];

    ASSERT_NEAR( predicted, expected, 1e-10 );
  }
}

void checkDirectionalDerivative( real64 const (&input)[4],
                                 real64 (& perturbedInput)[4],
                                 real64 const & val,
                                 real64 & perturbedVal,
                                 real64 const (&derivatives)[4],
                                 real64 (& perturbedDerivatives)[4],
                                 real64 const & perturb,
                                 real64 const & relTol,
                                 localIndex const direction,
                                 TableFunction::KernelWrapper kernelWrapper )
{
  LvArray::tensorOps::copy< 4 >( perturbedInput, input );
  real64 const dInput = perturb * ( input[direction] + perturb );
  perturbedInput[direction] += dInput;
  perturbedVal = kernelWrapper.compute( perturbedInput, perturbedDerivatives );

  geos::testing::checkRelativeError( derivatives[direction], (perturbedVal-val)/dInput, relTol, geos::testing::DEFAULT_ABS_TOL );
}

TEST( FunctionTests, 1DTable )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  // 1D table, various interpolation methods
  localIndex const Naxis = 4;
  localIndex const Ntest = 6;

  // Setup table
  array1d< array1d< real64 > > coordinates;
  coordinates.resize( 1 );
  coordinates[0].resize( Naxis );
  coordinates[0][0] = -1.0;
  coordinates[0][1] = 0.0;
  coordinates[0][2] = 2.0;
  coordinates[0][3] = 5.0;

  array1d< real64 > values( Naxis );
  values[0] = 1.0;
  values[1] = 3.0;
  values[2] = -5.0;
  values[3] = 7.0;

  TableFunction & table_a = dynamicCast< TableFunction & >( *functionManager->createChild( "TableFunction", "table_a" ) );
  table_a.setTableCoordinates( coordinates, { units::Dimensionless } );
  table_a.setTableValues( values, units::Dimensionless );
  table_a.reInitializeFunction();

  // Setup testing coordinates, expected values
  array1d< real64 > testCoordinates( Ntest );
  testCoordinates[0] = -1.1;
  testCoordinates[1] = -0.5;
  testCoordinates[2] = 0.2;
  testCoordinates[3] = 2.0;
  testCoordinates[4] = 4.0;
  testCoordinates[5] = 10.0;

  // Linear Interpolation
  array1d< real64 > testExpected( Ntest );
  testExpected[0] = 1.0;
  testExpected[1] = 2.0;
  testExpected[2] = 2.2;
  testExpected[3] = -5.0;
  testExpected[4] = 3.0;
  testExpected[5] = 7.0;
  table_a.setInterpolationMethod( TableFunction::InterpolationType::Linear );
  table_a.reInitializeFunction();
  evaluate1DFunction( table_a, testCoordinates, testExpected );

  // Upper
  testExpected[0] = 1.0;
  testExpected[1] = 3.0;
  testExpected[2] = -5.0;
  testExpected[3] = -5.0;
  testExpected[4] = 7.0;
  testExpected[5] = 7.0;
  table_a.setInterpolationMethod( TableFunction::InterpolationType::Upper );
  table_a.reInitializeFunction();
  evaluate1DFunction( table_a, testCoordinates, testExpected );

  // Lower
  testExpected[0] = 1.0;
  testExpected[1] = 1.0;
  testExpected[2] = 3.0;
  testExpected[3] = 3.0;
  testExpected[4] = -5.0;
  testExpected[5] = 7.0;
  table_a.setInterpolationMethod( TableFunction::InterpolationType::Lower );
  table_a.reInitializeFunction();
  evaluate1DFunction( table_a, testCoordinates, testExpected );

  // Nearest
  testExpected[0] = 1.0;
  testExpected[1] = 1.0;
  testExpected[2] = 3.0;
  testExpected[3] = -5.0;
  testExpected[4] = 7.0;
  testExpected[5] = 7.0;
  table_a.setInterpolationMethod( TableFunction::InterpolationType::Nearest );
  table_a.reInitializeFunction();
  evaluate1DFunction( table_a, testCoordinates, testExpected );

}



TEST( FunctionTests, 2DTable )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  // 2D table with linear interpolation
  // f(x, y) = 2*x - 3*y + 5
  localIndex const Ndim = 2;
  localIndex const Nx = 3;
  localIndex const Ny = 4;
  localIndex const Ntest = 20;
  string const inputName = "coordinates";

  // Setup table
  array1d< array1d< real64 > > coordinates;
  coordinates.resize( Ndim );
  coordinates[0].resize( Nx );
  coordinates[0][0] = -1.0;
  coordinates[0][1] = 0.0;
  coordinates[0][2] = 2.0;
  coordinates[1].resize( Ny );
  coordinates[1][0] = -1.0;
  coordinates[1][1] = 0.0;
  coordinates[1][2] = 1.0;
  coordinates[1][3] = 2.0;

  array1d< real64 > values( Nx * Ny );
  localIndex tablePosition = 0;
  for( localIndex jj=0; jj<Ny; ++jj )
  {
    for( localIndex ii=0; ii<Nx; ++ii )
    {
      real64 const x = coordinates[0][ii];
      real64 const y = coordinates[1][jj];
      values[tablePosition] = (2.0*x) - (3.0*y) + 5.0;
      ++tablePosition;
    }
  }

  // Set input names
  string_array inputVarNames( 1 );
  inputVarNames[0] = inputName;

  // Initialize the table
  TableFunction & table_b = dynamicCast< TableFunction & >( *functionManager->createChild( "TableFunction", "table_b" ) );
  table_b.setTableCoordinates( coordinates, { units::Dimensionless } );
  table_b.setTableValues( values, units::Dimensionless );
  table_b.setInterpolationMethod( TableFunction::InterpolationType::Linear );
  table_b.setInputVarNames( inputVarNames );
  table_b.reInitializeFunction();

  // Setup a group for testing the batch mode function evaluation
  conduit::Node node;
  dataRepository::Group testGroup( "testGroup", node );

  array2d< real64 > testCoordinates;
  testGroup.registerWrapper( inputName, &testCoordinates ).
    setSizedFromParent( 1 ).
    reference().resizeDimension< 1 >( Ndim );
  testGroup.resize( Ntest );

  // Build testing inputs/outputs
  array1d< real64 > expected( Ntest );
  array1d< real64 > output( Ntest );
  SortedArray< localIndex > set;

  // Build the set
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    set.insert( ii );
  }

  // Setup the random number generator
  std::default_random_engine generator;
  std::uniform_real_distribution< double > distribution( -0.99, 1.99 );

  // Test the function
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    for( localIndex jj=0; jj<Ndim; ++jj )
    {
      testCoordinates[ii][jj] = distribution( generator );
    }

    real64 const x = testCoordinates[ii][0];
    real64 const y = testCoordinates[ii][1];
    expected[ii] = (2.0*x) - (3.0*y) + 5.0;
    set.insert( ii );
  }

  // Evaluate the function in batch mode
  table_b.evaluate( testGroup, 0.0, set.toView(), output );

  // Compare results
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    ASSERT_NEAR( expected[ii], output[ii], 1e-10 );
  }
}


TEST( FunctionTests, 4DTable_multipleInputs )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  // 4D table with linear interpolation
  // f(x, y, z, t) = 2.0 + 3*x - 5*y + 7*z + 11*t
  localIndex const Ndim = 4;
  localIndex const Nx = 3;
  localIndex const Ny = 4;
  localIndex const Nz = 5;
  localIndex const Nt = 2;
  localIndex const Ntest = 20;
  localIndex const Ntimes = 5;
  string const coordinatesName = "coordinates";
  string const timeName = "time";

  // Setup table
  array1d< array1d< real64 > > coordinates;
  coordinates.resize( Ndim );
  coordinates[0].resize( Nx );
  coordinates[0][0] = -1.0;
  coordinates[0][1] = 0.0;
  coordinates[0][2] = 1.0;
  coordinates[1].resize( Ny );
  coordinates[1][0] = -1.0;
  coordinates[1][1] = 0.0;
  coordinates[1][2] = 0.5;
  coordinates[1][3] = 1.0;
  coordinates[2].resize( Nz );
  coordinates[2][0] = -1.0;
  coordinates[2][1] = -0.4;
  coordinates[2][2] = 0.3;
  coordinates[2][3] = 0.5;
  coordinates[2][4] = 1.0;
  coordinates[3].resize( Nt );
  coordinates[3][0] = -1.0;
  coordinates[3][1] = 1.0;

  array1d< real64 > values( Nx * Ny * Nz * Nt );
  localIndex tablePosition = 0;
  for( localIndex mm=0; mm<Nt; ++mm )
  {
    for( localIndex kk=0; kk<Nz; ++kk )
    {
      for( localIndex jj=0; jj<Ny; ++jj )
      {
        for( localIndex ii=0; ii<Nx; ++ii )
        {
          real64 const x = coordinates[0][ii];
          real64 const y = coordinates[1][jj];
          real64 const z = coordinates[2][kk];
          real64 const t = coordinates[3][mm];
          values[tablePosition] = 2.0 + (3*x) - (5*y) + (7*z) + (11*t);
          ++tablePosition;
        }
      }
    }
  }

  // Set variable names
  string_array inputVarNames( 2 );
  inputVarNames[0] = coordinatesName;
  inputVarNames[1] = timeName;

  // Initialize the table
  TableFunction & table_c = dynamicCast< TableFunction & >( *functionManager->createChild( "TableFunction", "table_c" ) );
  table_c.setTableCoordinates( coordinates, { units::Dimensionless } );
  table_c.setTableValues( values, units::Dimensionless );
  table_c.setInterpolationMethod( TableFunction::InterpolationType::Linear );
  table_c.setInputVarNames( inputVarNames );
  table_c.reInitializeFunction();

  // Setup a group for testing the batch mode function evaluation
  conduit::Node node;
  dataRepository::Group testGroup( "testGroup", node );

  array2d< real64 > testCoordinates;
  testGroup.registerWrapper( coordinatesName, &testCoordinates ).
    setSizedFromParent( 1 ).
    reference().resizeDimension< 1 >( Ndim - 1 );
  testGroup.resize( Ntest );

  // Build testing inputs/outputs
  array1d< real64 > expected( Ntest );
  array1d< real64 > output( Ntest );
  SortedArray< localIndex > set;

  // Fill out the set
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    set.insert( ii );
  }

  // Setup a random number generator
  std::default_random_engine generator;
  std::uniform_real_distribution< double > distribution( -0.99, 0.99 );

  // Build the inputs
  for( localIndex ii=0; ii<Ntimes; ++ii )
  {
    real64 const t = distribution( generator );

    for( localIndex jj=0; jj<Ntest; ++jj )
    {
      for( localIndex kk=0; kk<Ndim-1; ++kk )
      {
        testCoordinates[jj][kk] = distribution( generator );
      }

      real64 const x = testCoordinates[jj][0];
      real64 const y = testCoordinates[jj][1];
      real64 const z = testCoordinates[jj][2];
      expected[jj] = 2.0 + (3*x) - (5*y) + (7*z) + (11*t);
    }

    // Evaluate the function in batch mode
    table_c.evaluate( testGroup, t, set.toView(), output );

    // Compare results
    for( localIndex jj=0; jj<Ntest; ++jj )
    {
      ASSERT_NEAR( expected[jj], output[jj], 1e-10 );
    }
  }
}

TEST( FunctionTests, 4DTable_derivatives )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  // 4D table with linear interpolation
  // f(x, y, z, t) = 2.0 + 3*x - 5*y*y + 7*z*z*z + 11*t*t*t*t
  localIndex const Ndim = 4;
  localIndex const Nx = 3;
  localIndex const Ny = 4;
  localIndex const Nz = 5;
  localIndex const Nt = 4;
  string const coordinatesName = "coordinates";
  string const timeName = "time";

  // Setup table
  array1d< array1d< real64 > > coordinates;
  coordinates.resize( Ndim );
  coordinates[0].resize( Nx );
  coordinates[0][0] = -1.0;
  coordinates[0][1] = 0.0;
  coordinates[0][2] = 1.0;
  coordinates[1].resize( Ny );
  coordinates[1][0] = -1.0;
  coordinates[1][1] = 0.0;
  coordinates[1][2] = 0.5;
  coordinates[1][3] = 1.0;
  coordinates[2].resize( Nz );
  coordinates[2][0] = -1.0;
  coordinates[2][1] = -0.4;
  coordinates[2][2] = 0.3;
  coordinates[2][3] = 0.5;
  coordinates[2][4] = 1.0;
  coordinates[3].resize( Nt );
  coordinates[3][0] = -1.0;
  coordinates[3][1] = 0.34;
  coordinates[3][2] = 0.5;
  coordinates[3][3] = 1.0;

  array1d< real64 > values( Nx * Ny * Nz * Nt );
  for( localIndex mm=0, tablePosition=0; mm<Nt; ++mm )
  {
    for( localIndex kk=0; kk<Nz; ++kk )
    {
      for( localIndex jj=0; jj<Ny; ++jj )
      {
        for( localIndex ii=0; ii<Nx; ++ii, ++tablePosition )
        {
          real64 const x = coordinates[0][ii];
          real64 const y = coordinates[1][jj];
          real64 const z = coordinates[2][kk];
          real64 const t = coordinates[3][mm];

          values[tablePosition] = 2.0 + (3*x) - (5*y*y) + (7*z*z*z) + (11*t*t*t*t);
        }
      }
    }
  }

  // Set variable names
  string_array inputVarNames( 2 );
  inputVarNames[0] = coordinatesName;
  inputVarNames[1] = timeName;

  // Initialize the table
  TableFunction & table_d = dynamicCast< TableFunction & >( *functionManager->createChild( "TableFunction", "table_d" ) );
  table_d.setTableCoordinates( coordinates, { units::Dimensionless } );
  table_d.setTableValues( values, units::Dimensionless );
  table_d.setInterpolationMethod( TableFunction::InterpolationType::Linear );
  table_d.setInputVarNames( inputVarNames );
  table_d.reInitializeFunction();

  real64 const relTol = 1e-4;
  real64 const eps = std::numeric_limits< real64 >::epsilon();
  real64 const perturb = std::sqrt( eps );

  real64 const start = coordinates[0][0]-0.1; // start outside the table, to check derivatives there too
  real64 const end = coordinates[0][Nx-1]+0.1; // end outside the table
  localIndex const nSamples = 7; // try not to fall on table coordinate, otherwise the finite-difference approximation won't work
  real64 const delta = (end-start)/nSamples;

  real64 val = 0.0;
  real64 perturbedVal = 0.0;
  real64 input[4] = { start, start, start, start };
  real64 perturbedInput[4]{};
  real64 derivatives[4]{};
  real64 perturbedDerivatives[4]{};

  TableFunction::KernelWrapper kernelWrapper = table_d.createKernelWrapper();
  for( localIndex mm=0; mm<nSamples; ++mm, input[3] += delta )
  {
    input[2] = start;
    for( localIndex kk=0; kk<nSamples; ++kk, input[2] += delta )
    {
      input[1] = start;
      for( localIndex jj=0; jj<nSamples; ++jj, input[1] += delta )
      {
        input[0] = start;
        for( localIndex ii=0; ii<nSamples; ++ii, input[0] += delta )
        {
          // evaluate once to get the analytical derivatives
          val = kernelWrapper.compute( input, derivatives );

          // check derivatives

          for( localIndex direction = 0; direction < Ndim; ++direction )
          {
            checkDirectionalDerivative( input, perturbedInput,
                                        val, perturbedVal,
                                        derivatives, perturbedDerivatives,
                                        perturb, relTol, direction, kernelWrapper );
          }
        }
      }
    }
  }
}

#ifdef GEOS_USE_MATHPRESSO

TEST( FunctionTests, 4DTable_symbolic )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  // Symbolic function with four inputs
  string const expression = "1.0+(2.0*a)-(3.0*b*b)+(5.0*c*c*c)-(7.0*d*d*d*d)";
  localIndex const Ntest = 20;

  // Set variable names
  string_array inputVarNames( 4 );
  string const nameA = "a";
  string const nameB = "b";
  string const nameC = "c";
  string const nameD = "d";
  inputVarNames[0] = nameA;
  inputVarNames[1] = nameB;
  inputVarNames[2] = nameC;
  inputVarNames[3] = nameD;

  // Initialize the table
  SymbolicFunction & table_e = dynamicCast< SymbolicFunction & >( *functionManager->createChild( "SymbolicFunction", "table_e" ) );
  table_e.setSymbolicExpression( expression );
  table_e.setInputVarNames( inputVarNames );
  table_e.setSymbolicVariableNames( inputVarNames );
  table_e.initializeFunction();

  // Setup a group for testing the batch mode function evaluation
  conduit::Node node;
  dataRepository::Group testGroup( "testGroup", node );

  array1d< real64 > inputA;
  array1d< real64 > inputB;
  array1d< real64 > inputC;
  array1d< real64 > inputD;
  testGroup.registerWrapper( nameA, &inputA ).setSizedFromParent( 1 );
  testGroup.registerWrapper( nameB, &inputB ).setSizedFromParent( 1 );
  testGroup.registerWrapper( nameC, &inputC ).setSizedFromParent( 1 );
  testGroup.registerWrapper( nameD, &inputD ).setSizedFromParent( 1 );
  testGroup.resize( Ntest );

  // Build testing inputs/outputs
  array1d< real64 > expected( Ntest );
  array1d< real64 > output( Ntest );
  SortedArray< localIndex > set;

  // Fill out the set
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    set.insert( ii );
  }

  // Setup a random number generator
  std::default_random_engine generator;
  std::uniform_real_distribution< double > distribution( -1.0, 1.0 );

  // Build the inputs
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    real64 const a = distribution( generator );
    real64 const b = distribution( generator );
    real64 const c = distribution( generator );
    real64 const d = distribution( generator );
    inputA[ii] = a;
    inputB[ii] = b;
    inputC[ii] = c;
    inputD[ii] = d;

    expected[ii] = 1.0+(2.0*a)-(3.0*b*b)+(5.0*c*c*c)-(7.0*d*d*d*d);
  }

  // Evaluate the function in batch mode
  table_e.evaluate( testGroup, 0.0, set.toView(), output );

  // Compare results
  for( localIndex jj=0; jj<Ntest; ++jj )
  {
    ASSERT_NEAR( expected[jj], output[jj], 1e-10 );
  }
}

#endif

template< integer NUM_DIMS, integer NUM_OPS >
void testMutivariableFunction( MultivariableTableFunction & function,
                               arrayView1d< real64 const > const & inputs,
                               arrayView1d< real64 const > const & expectedValues,
                               arrayView1d< real64 const > const & expectedDerivatives,
                               real64 valuesTolerance = 1e-10,
                               real64 derivativesTolerance = 1e-10 )
{
  localIndex const numElems = inputs.size() / NUM_DIMS;

  ASSERT_EQ( numElems * NUM_DIMS, inputs.size());
  ASSERT_EQ( numElems * NUM_OPS, expectedValues.size());

  array1d< real64 > evaluatedValues( numElems * NUM_OPS );
  array2d< real64 > evaluatedDerivatives( numElems * NUM_OPS, NUM_DIMS );
  arrayView1d< real64 > evaluatedValuesView = evaluatedValues.toView();
  arrayView2d< real64 > evaluatedDerivativesView = evaluatedDerivatives.toView();


  MultivariableTableFunctionStaticKernel< NUM_DIMS, NUM_OPS > kernel( function.getAxisMinimums(),
                                                                      function.getAxisMaximums(),
                                                                      function.getAxisPoints(),
                                                                      function.getAxisSteps(),
                                                                      function.getAxisStepInvs(),
                                                                      function.getAxisHypercubeMults(),
                                                                      function.getHypercubeData()
                                                                      );
  // Test values evaluation first
  forAll< geos::parallelDevicePolicy< > >( numElems, [=] GEOS_HOST_DEVICE
                                             ( localIndex const elemIndex )
  {
    kernel.compute( &inputs[elemIndex * NUM_DIMS], &evaluatedValuesView[elemIndex * NUM_OPS] );
  } );

  forAll< serialPolicy >( numElems * NUM_OPS, [=] ( localIndex const elemOpIndex )
  {
    ASSERT_NEAR( expectedValues[elemOpIndex], evaluatedValuesView[elemOpIndex], valuesTolerance );
  } );

  // And now - both values and derivatives
  forAll< geos::parallelDevicePolicy< > >( numElems, [=] GEOS_HOST_DEVICE
                                             ( localIndex const elemIndex )
  {
    // use local 2D array for the kernel
    real64 derivatives[NUM_OPS][NUM_DIMS];

    kernel.compute( &inputs[elemIndex * NUM_DIMS], &evaluatedValuesView[elemIndex * NUM_OPS], derivatives );

    // now copy results to the view
    for( auto i = 0; i < NUM_OPS; i++ )
      for( auto j = 0; j < NUM_DIMS; j++ )
        evaluatedDerivativesView[elemIndex * NUM_OPS + i][j] = derivatives[i][j];
  } );

  // Perform checks.
  forAll< serialPolicy >( numElems * NUM_OPS, [=] ( localIndex const elemOpIndex )
  {
    ASSERT_NEAR( expectedValues[elemOpIndex], evaluatedValuesView[elemOpIndex], valuesTolerance );
  } );

  // Perform checks.
  forAll< serialPolicy >( numElems * NUM_OPS, [=] ( localIndex const elemOpIndex )
  {
    for( auto j = 0; j < NUM_DIMS; j++ )
      ASSERT_NEAR( expectedDerivatives[elemOpIndex * NUM_DIMS + j], evaluatedDerivativesView[elemOpIndex][j], derivativesTolerance );
  } );
}

TEST( FunctionTests, 1DMultivariableTable )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  // 1D table
  localIndex constexpr nDims = 1;
  localIndex constexpr nOps = 1;
  localIndex const nTest = 7;

  // Setup table
  array1d< real64 > axisMins;
  array1d< real64 > axisMaxs;
  integer_array axisPoints;
  axisMins.resize( nDims );
  axisMaxs.resize( nDims );
  axisPoints.resize( nDims );
  axisMins[0] = 1;
  axisMaxs[0] = 100;
  axisPoints[0] = 100;

  array1d< real64 > values( axisPoints[0] );
  // set this up so that values are equal to coordinate: f(x) = x
  for( auto i = 0; i < axisPoints[0]; i++ )
    values[i] = i + 1;

  MultivariableTableFunction & table_f = dynamicCast< MultivariableTableFunction & >( *functionManager->createChild( "MultivariableTableFunction", "table_f" ) );
  table_f.setTableCoordinates( nDims, nOps, axisMins, axisMaxs, axisPoints );
  table_f.setTableValues( values );
  table_f.initializeFunction();


  // Setup testing coordinates, expected values
  array1d< real64 > testCoordinates( nTest );
  testCoordinates[0] = -1.1;
  testCoordinates[1] = 0.0;
  testCoordinates[2] = 1.0;
  testCoordinates[3] = 2.8;
  testCoordinates[4] = 99.3;
  testCoordinates[5] = 100.0;
  testCoordinates[6] = 11234.534;

  array1d< real64 > testExpectedValues( testCoordinates );
  array1d< real64 > testExpectedDerivatives( nTest * nDims * nOps );
  for( auto i = 0; i < nTest * nDims * nOps; i++ )
    testExpectedDerivatives[i] = 1;


  testMutivariableFunction< nDims, nOps >( table_f, testCoordinates, testExpectedValues, testExpectedDerivatives );
}

real64 operator1 ( real64 const x, real64 const y ) { return 2 * x + 3 * y * y; }
real64 dOperator1_dx ( real64 const x, real64 const y ) { GEOS_UNUSED_VAR( x, y ); return 2; }
real64 dOperator1_dy ( real64 const x, real64 const y ) { GEOS_UNUSED_VAR( x ); return 6 * y; }
real64 operator2 ( real64 const x, real64 const y ) { return 2 * x * x + 3 / (y + 1); }
real64 dOperator2_dx ( real64 const x, real64 const y ) { GEOS_UNUSED_VAR( y ); return 4 * x; }
real64 dOperator2_dy ( real64 const x, real64 const y ) { GEOS_UNUSED_VAR( x ); return -3 / ((y + 1) * (y + 1)); }
real64 operator3 ( real64 const x, real64 const y ) { GEOS_UNUSED_VAR( y ); return 2 * x + 3; }
real64 dOperator3_dx ( real64 const x, real64 const y ) { GEOS_UNUSED_VAR( x, y ); return 2; }
real64 dOperator3_dy ( real64 const x, real64 const y ) { GEOS_UNUSED_VAR( x, y ); return 0; }

TEST( FunctionTests, 2DMultivariableTable )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();


  // 1D table
  localIndex constexpr nDims = 2;
  localIndex constexpr nOps = 3;
  localIndex const nTest = 3;

  // Setup table
  array1d< real64 > axisMins( nDims );
  array1d< real64 > axisMaxs( nDims );
  integer_array axisPoints( nDims );


  axisMins[0] = 1;
  axisMins[1] = 0;
  axisMaxs[0] = 2;
  axisMaxs[1] = 1;
  axisPoints[0] = 1000;
  axisPoints[1] = 1100;

  array1d< real64 > axisSteps( nDims );
  for( auto i = 0; i < nDims; i++ )
    axisSteps[i] = (axisMaxs[i] - axisMins[i]) /  axisPoints[i];

  array1d< real64 > values( axisPoints[0] * axisPoints[1] * nOps );
  for( auto i = 0; i < axisPoints[0]; i++ )
    for( auto j = 0; j < axisPoints[1]; j++ )
    {
      auto x = axisMins[0] + i * axisSteps[0];
      auto y = axisMins[1] + j * axisSteps[1];

      values[( i * axisPoints[1] + j ) * nOps] = operator1( x, y );
      values[( i * axisPoints[1] + j ) * nOps + 1] = operator2( x, y );
      values[( i * axisPoints[1] + j ) * nOps + 2] = operator3( x, y );
    }

  MultivariableTableFunction & table_g = dynamicCast< MultivariableTableFunction & >( *functionManager->createChild( "MultivariableTableFunction", "table_f" ) );
  table_g.setTableCoordinates( nDims, nOps, axisMins, axisMaxs, axisPoints );
  table_g.setTableValues( values );
  table_g.initializeFunction();


  // Setup testing coordinates, expected values
  array1d< real64 > testCoordinates( nTest * nDims );
  testCoordinates[0] = 1.2334;
  testCoordinates[1] = 0.1232;
  testCoordinates[2] = 1.7342;
  testCoordinates[3] = 0.2454;
  testCoordinates[4] = 2.0;
  testCoordinates[5] = 0.7745;

  array1d< real64 > testExpectedValues( nTest * nOps );
  array1d< real64 > testExpectedDerivatives( nTest * nOps * nDims );
  for( auto i = 0; i < nTest; i++ )
  {
    testExpectedValues[i * nOps] = operator1( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedValues[i * nOps + 1] = operator2( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedValues[i * nOps + 2] = operator3( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedDerivatives[i * nOps * nDims] = dOperator1_dx( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedDerivatives[i * nOps * nDims + 1] = dOperator1_dy( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedDerivatives[i * nOps * nDims + 2] = dOperator2_dx( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedDerivatives[i * nOps * nDims + 3] = dOperator2_dy( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedDerivatives[i * nOps * nDims + 4] = dOperator3_dx( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedDerivatives[i * nOps * nDims + 5] = dOperator3_dy( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
  }


  testMutivariableFunction< nDims, nOps >( table_g, testCoordinates, testExpectedValues, testExpectedDerivatives, 1e-2, 2e-2 );
}

TEST( FunctionTests, MultivariableTableFromFile )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  localIndex const nTest = 1;
  localIndex constexpr nDims = 4;
  localIndex constexpr nOps = 34;

  // Setup table

  MultivariableTableFunction & table_h = dynamicCast< MultivariableTableFunction & >( *functionManager->createChild( "MultivariableTableFunction", "table_f" ) );

  writeTableToFile( "tableData.txt", multivariableTableFileContent );
  table_h.initializeFunctionFromFile ( "tableData.txt" );
  removeFile( "tableData.txt" );


  // Setup testing coordinates, expected values
  array1d< real64 > testCoordinates( nTest * 4 );
  testCoordinates[0] = 345;
  testCoordinates[1] = 0.1232;
  testCoordinates[2] = 0.2454;
  testCoordinates[3] = 0.4745;

  real64 vals[nTest *
              nOps] =
  { 2.425619913180601, 8.036906715392849, 30.043211252986385, 3.135411493315748, 0.07018930898607319, 6.627254956030642, 27.89619715347374, 9.51145988640026e-13, 29.178033463062807,
    9.962409391432278e-12, 0.12583684837156695, 4.708882886399814e-13, 0.674544852574673, 0.16858415006532698, 0.00020186000309950752, 0.0028632490461658346, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -1.1408215389715537, 0.0, 0.5704107694857768, 0.0, 0.0, 0.0, 1340.4622053037895, 292.252107461518, 0.0, 0.0, 0.8431 };

  array1d< real64 > testExpectedValues ( nTest * nOps );
  for( auto i = 0; i < nTest * nOps; i++ )
  {
    testExpectedValues[i] = vals[i];
  }

  real64 ders[nTest *
              nOps* nDims] =
  {0.0001940148221947931, 19.68847332124942, 0.4062830261664956, 1.1366559830488714, 4.8455276643095086e-05, -2.9592626628654695, 32.75023111406469, 3.0346605840599614, 0.0001393383426014372,
   -20.98597627326845, -15.133478247625717, 63.31551370493199, 3.135303638867035e-07, -19.983502188118216, -19.983502188118216, -19.98350218811822, 7.964896387703173e-06, -0.023281125656838278,
   -0.07157042919964876, -0.11429929487246732, 0.0001581574364548, -6.598332756805898, 20.21093196228746, -3.9775257787170633, 0.00018582807579339933, -30.549324995402156, -25.34377155909867,
   53.767847848009815, 0.0, -3.965423000001099e-13, 6.4741599999974e-14, 9.296671999996806e-14, 0.0012618474420348594, 86.03146322988385, -29.81083896358708, -35.34918506365271,
   2.322189453563556e-16, 6.806791366159914e-12, -6.324028712102855e-12, -1.763423276482453e-11, 4.194925057674508e-06, 0.28426943301429053, -0.1490178794631717, -0.054137613393434, 0.0,
   6.034577000000902e-13, -4.607584000001362e-13, -6.616332800002337e-13, 2.9147682999364393e-06, -0.36834344143989595, 0.9765255961965139, 0.9244472245352091, -2.830458299936535e-06,
   1.368377841439896, 0.023508803803486127, 0.07558717546479085, 2.0179058713752443e-10, 0.00020652903147683487, -0.00020665974970602544, -0.0002145141199359456, 2.7680288536019274e-07,
   0.003671016958495129, -0.002800695264266917, -0.004025376506406397, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   0.0, 0.0, 0.0, 0.0, 5.3196256172338705, 5.047052119868672, 1.8075747560203723, 0.0, 0.0, 0.0, 0.0, 0.0, -2.6598128086169353, -2.523526059934336, -0.9037873780101862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0013400012448758296, -551.2530669025564, 82.79691651792007, 118.86424502561012, 0.028253297318398797, 374.5299869280559, -285.9650933760846, -410.6360788993451,
   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };


  array1d< real64 > testExpectedDerivatives( nTest * nOps * nDims );

  for( auto i = 0; i < nTest * nOps * nDims; i++ )
  {
    testExpectedDerivatives[i] = ders[i];
  }

  testMutivariableFunction< nDims, nOps >( table_h, testCoordinates, testExpectedValues, testExpectedDerivatives );
}


// The `ENUM_STRING` implementation relies on consistency between the order of the `enum`,
// and the order of the `string` array provided. Since this consistency is not enforced, it can be corrupted anytime.
// This unit test aims at preventing from this implicit relationship to bring a bug.
TEST( TableFunctionEnums, InterpolationType )
{
  using EnumType = TableFunction::InterpolationType;

  ASSERT_EQ( "linear", toString( EnumType::Linear ) );
  ASSERT_EQ( "nearest", toString( EnumType::Nearest ) );
  ASSERT_EQ( "upper", toString( EnumType::Upper ) );
  ASSERT_EQ( "lower", toString( EnumType::Lower ) );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  // geos::GeosxState state( geos::basicSetup( argc, argv ) );
  conduit::Node conduitNode;
  dataRepository::Group rootNode( "root", conduitNode );


  FunctionManager functionManager( "FunctionManager", &rootNode );

  int const result = RUN_ALL_TESTS();

  // geos::basicCleanup();

  return result;
}
