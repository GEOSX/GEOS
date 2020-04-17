#include "VTKTypes.hpp"

#include "common/DataTypes.hpp"

#include <vtkType.h>
#include <vtkCellType.h>

#include <unordered_map>

namespace geosx
{
namespace vtk
{

static std::unordered_map< std::type_index, int > geosxToVTKTypeMap =
{
  {std::type_index( typeid( integer ) ), VTK_INT},
  {std::type_index( typeid( localIndex ) ), VTK_INT},
  {std::type_index( typeid( globalIndex ) ), VTK_LONG},
  {std::type_index( typeid( real32 ) ), VTK_FLOAT},
  {std::type_index( typeid( real64 ) ), VTK_DOUBLE}
};

static std::unordered_map< string, int > geosxToVTKCellTypeMap =
{
  { "C3D4", VTK_HEXAHEDRON },
  { "C3D5", VTK_PYRAMID },
  { "C3D6", VTK_WEDGE },
  { "C3D8", VTK_HEXAHEDRON }
};

int GetVTKType( std::type_index const & type )
{
  int VTKType = -1;
  try
  {
    VTKType = geosxToVTKTypeMap.at( type );
  }
  catch( const std::out_of_range& outOfRange )
  {
  GEOSX_ERROR( "Type " << type.name() << " not recognized for a VTK output" );
  }
  return VTKType;
}

int GetVTKCellType( string const & type )
{
  int VTKCellType = -1;
  try
  {
    VTKCellType = geosxToVTKCellTypeMap.at( type );
  }
  catch( const std::out_of_range& outOfRange )
  {
    GEOSX_ERROR( "Cell type " << type << " not recognized for a VTK output" );
  }
  return VTKCellType;
}
}
}
