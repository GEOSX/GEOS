

#include "physicsSolvers/solidMechanics/kernels/SmallStrainResidual_impl.hpp"
#include "policies.hpp"



#define INSTANTIATION( NAME ) \
  template class NAME< CellElementSubRegion, ElasticIsotropic, H1_Hexahedron_Lagrange1_GaussLegendre2 >; \
  template real64 NAME< CellElementSubRegion, ElasticIsotropic, H1_Hexahedron_Lagrange1_GaussLegendre2 >::kernelLaunch< NAME ## Policy, \
                                                                                                                        NAME< CellElementSubRegion, ElasticIsotropic, \
                                                                                                                              H1_Hexahedron_Lagrange1_GaussLegendre2 > > \
    ( localIndex const, \
    NAME< CellElementSubRegion, ElasticIsotropic, H1_Hexahedron_Lagrange1_GaussLegendre2 > const & ); \


namespace geos
{
using namespace constitutive;
using namespace finiteElement;
namespace solidMechanicsLagrangianFEMKernels
{
INSTANTIATION( SmallStrainResidual )
}
}

#undef INSTANTIATION
