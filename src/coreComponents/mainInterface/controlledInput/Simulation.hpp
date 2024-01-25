#ifndef GEOS_SIMULATION_HPP
#define GEOS_SIMULATION_HPP

#include "Solvers.hpp"

#include "common/DataTypes.hpp"

#include <pugixml.hpp>
#include <yaml-cpp/yaml.h>

namespace geos::input
{

class Simulation
{
public:

  void setBegin( string const & begin );

  void setEnd( string const & end );

  void setSolver( std::vector< std::shared_ptr< solvers::Solver > > const & solvers );

  void fillProblemXmlNode( pugi::xml_node & problemNode ) const;

private:

  void fillEventsXmlNode( pugi::xml_node & eventsNode ) const;

  string m_begin;
  string m_end;
  std::vector< std::shared_ptr< solvers::Solver > > m_solvers;
};

void operator>>( const YAML::Node & node,
                 Simulation & simulation );

}

#endif //GEOS_SIMULATION_HPP
