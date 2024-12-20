
#include "common/GeosxMacros.hpp"

// TPL includes
#include <gtest/gtest.h>


namespace geos
{
namespace denseLinearAlgebra
{
namespace testing
{


#if defined(GEOS_DEVICE_COMPILE)
#define PORTABLE_EXPECT_EQ( L, R ) GEOS_ERROR_IF_NE( L, R )
#define PORTABLE_EXPECT_NEAR( L, R, EPSILON ) LVARRAY_ERROR_IF_GE_MSG( LvArray::math::abs( ( L ) -( R ) ), EPSILON, \
                                                                       STRINGIZE( L ) " = " << ( L ) << "\n" << STRINGIZE( R ) " = " << ( R ) );
#define PORTABLE_EXPECT_TRUE( value ) GEOS_ERROR_IF( !value, "should be true" )
#define PORTABLE_EXPECT_FALSE( value ) GEOS_ERROR_IF( value, "should be false" )
#else
#define PORTABLE_EXPECT_EQ( L, R ) EXPECT_EQ( L, R )
#define PORTABLE_EXPECT_NEAR( L, R, EPSILON ) EXPECT_LE( LvArray::math::abs( ( L ) -( R ) ), EPSILON ) << \
    STRINGIZE( L ) " = " << ( L ) << "\n" << STRINGIZE( R ) " = " << ( R );
#define PORTABLE_EXPECT_TRUE( value ) EXPECT_TRUE( value )
#define PORTABLE_EXPECT_FALSE( value ) EXPECT_FALSE( value )
#endif

} //namespace testing

} // namespace denseLinearAlgebra

} // namespace geos
