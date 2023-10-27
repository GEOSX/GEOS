#include "ControlledInput.hpp"

#include "dataRepository/xmlWrapper.hpp"

#include "codingUtilities/StringUtilities.hpp"

#include <yaml-cpp/yaml.h>

#include <map>
#include <vector>

namespace geos::input {

using namespace pugi;

real64 convertTime(string const & time)
{
  std::vector< string > tokens = stringutilities::tokenize< std::vector >( time, " " );
  std::size_t const numTokens = tokens.size();
  GEOS_ASSERT( numTokens == 1 || numTokens == 2 );

  real64 const t = std::stod( tokens.front() ); // TODO check cast.

  real64 scaling = 1.;
  if( numTokens == 2 )
  {
    GEOS_ASSERT_EQ( tokens.size(), 2 );
    int constexpr sec = 1;
    int constexpr min = 60;
    int constexpr hour = 60 * min;
    int constexpr day = 24 * hour;
    int constexpr week = 7 * day;
    real64 constexpr year = 365.25 * day;
    real64 constexpr month = year / 12.;

    std::map< string, real64 > const conv{  // TODO To real64 because of the 365.25
      { "s",       sec },
      { "sec",     sec },
      { "second",  sec },
      { "second(s)",  sec },
      { "seconds", sec },
      { "min",     min },
      { "minute",  min },
      { "minute(s)",  min },
      { "minutes", min },
      { "h",       hour },
      { "hour",    hour },
      { "hour(s)",    hour },
      { "hours",   hour },
      { "d",       day },
      { "day",     day },
      { "day(s)",     day },
      { "days",    day },
      { "w",       week },
      { "week",    week },
      { "week(s)",    week },
      { "weeks",   week },
      { "m",       month },
      { "month",   month },
      { "month(s)",   month },
      { "months",  month },
      { "y",       year },
      { "year",    year },
      { "year(s)",    year },
      { "years",   year },
    };
    scaling = conv.at( tokens.back() );// TODO check if value is found
  }
  return t * scaling;
}

class Simulation
{
public:

  void setBegin( string const & begin )
  {
    m_begin = begin;
  }

  void setEnd( string const & end )
  {
    m_end = end;
  }

  void fillEventsXmlNode( xml_node & eventsNode ) const
  {
    eventsNode.append_attribute( "minTime" ) = convertTime( m_begin );
    eventsNode.append_attribute( "maxTime" ) = convertTime( m_end );
  }

private:
  string m_begin;
  string m_end;
};

class Event
{
public:
  explicit Event( string const & target )
    : m_target( target )
  { }

  virtual ~Event() = default;

  virtual void fillEventsXmlNode( xml_node & eventsNode ) const = 0;
protected:
  string m_target;
};

class PeriodicEvent : public Event
{
public:
  PeriodicEvent( string const & target, string const & every )
    : Event( target ),
      m_every( every )
  { }

  void fillEventsXmlNode( xml_node & eventsNode ) const override
  {
    xml_node periodicEvent = eventsNode.append_child( "PeriodicEvent" );
    periodicEvent.append_attribute( "timeFrequency" ) = convertTime( m_every );
    periodicEvent.append_attribute( "target" ) = m_target.c_str();
  }

private:
  string m_every;
};

class SoloEvent : public Event
{
public:
  SoloEvent( string const & target, string const & at )
    : Event( target ),
      m_at( at )
  { }

  void fillEventsXmlNode( xml_node & eventsNode ) const override
  {
    xml_node soloEvent = eventsNode.append_child( "SoloEvent" );
    soloEvent.append_attribute( "targetTime" ) = convertTime( m_at );
    soloEvent.append_attribute( "target" ) = m_target.c_str();
  }

private:
  string m_at;
};

class Output
{
public:
  Output( string name,
          string const & every,
          std::vector< string > const & at )
    : m_name( name ),
      m_every( every ),
      m_at( at )
  { }

  virtual ~Output() = default;

  virtual void fillOutputsXmlNode( xml_node & outputsNode ) const = 0;

  std::vector< std::shared_ptr< Event > > getEvents() const
  {
    std::vector< std::shared_ptr< Event > > result;

    if( !m_every.empty() )
    {
      result.emplace_back( std::make_shared< PeriodicEvent >( "/Outputs/" + m_name, m_every ) );
    }

    for( string const & at: m_at )
    {
      if( !at.empty() )
      {
        result.emplace_back( std::make_shared< SoloEvent >( "/Outputs/" + m_name, at ) );
      }
    }

    return result;
  };

protected:
  string m_name;

private:
  string m_every;
  std::vector< string > m_at;
};

class VtkOutput: public Output
{
public:
  VtkOutput( int counter,
             string const & every,
             std::vector< string > const & at,
             std::vector< string > const & fields,
             bool writeGhostCells,
             string fileRoot )
    : Output( "__vtk-" + std::to_string( counter ), every, at ),
      m_fields( fields ),
      m_writeGhostCells( writeGhostCells ),
      m_fileRoot( fileRoot )
  { }

  void fillOutputsXmlNode( xml_node & outputsNode ) const override
  {
    xml_node vtkOutput = outputsNode.append_child( "VTK" );
    vtkOutput.append_attribute( "name" ) = m_name.c_str();
    vtkOutput.append_attribute( "writeGhostCells" ) = m_writeGhostCells ? "1" : "0";
    vtkOutput.append_attribute( "fieldNames" ) = ( "{ " + stringutilities::join( m_fields, ", " ) + " }" ).c_str();
    if( !m_fileRoot.empty() )
    {
      vtkOutput.append_attribute( "plotFileRoot" ) = m_fileRoot.c_str();
    }
  }

private:
  std::vector< string > m_fields;
  bool m_writeGhostCells;
  string m_fileRoot;
};


class RestartOutput: public Output
{
public:
  RestartOutput( int counter,
                 string const & every,
                 std::vector< string > const & at )
    : Output( "__restart-" + std::to_string( counter ), every, at )
  { }

  void fillOutputsXmlNode( xml_node & outputsNode ) const override
  {
    xml_node vtkOutput = outputsNode.append_child( "Restart" );
    vtkOutput.append_attribute( "name" ) = m_name.c_str();
  }
};

class Deck
{
public:

  void setSimulation( Simulation const & simulation )
  {
    m_simulation = simulation;
  }

  void setOutputs( std::vector< std::shared_ptr< Output > > const & outputs )
  {
    m_outputs = outputs;
  }

  void fillProblemXmlNode( xml_node & problemNode ) const
  {
    xml_node xmlOutputs = problemNode.append_child( "Outputs" );
    xml_node xmlEvents = problemNode.append_child( "Events" );

    m_simulation.fillEventsXmlNode( xmlEvents );

    for( std::shared_ptr< Output > output: m_outputs )
    {
      output->fillOutputsXmlNode( xmlOutputs );
      for( std::shared_ptr< Event > event: output->getEvents() )
      {
        event->fillEventsXmlNode( xmlEvents );
      }
    }

    // Add name to all the events.
    int iEvent = 0;
    for( auto it = xmlEvents.children().begin(); it != xmlEvents.children().end(); ++it, ++iEvent )
    {
      it->append_attribute( "name" ) = ( "__event-" + std::to_string( iEvent ) ).c_str();
    }
  }

private:
  Simulation m_simulation;
  std::vector< std::shared_ptr< Output > > m_outputs;
};

void operator>>( const YAML::Node & node,
                 Simulation & simulation )
{
  simulation.setBegin( node["begin"].as< string >() );
  simulation.setEnd( node["end"].as< string >() );
}

void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< Output > > & outputs )
{
  GEOS_ASSERT( node.IsSequence() );

  for( std::size_t i = 0; i < node.size(); i++ )
  {
    auto output = node[i];
    for (auto const & kv: output)
    {
      string const outputType = kv.first.as< string >();
      auto const & subNode = kv.second;
      if( outputType == "vtk" )
      {
        auto sp = std::make_shared< VtkOutput >( i,
                                                 subNode["every"].as< string >( "" ),
                                                 subNode["at"].as< std::vector< string > >( std::vector< string >() ),
                                                 subNode["fields"].as< std::vector< string > >(),
                                                 subNode["write_ghost_cells"].as< bool >(),
                                                 subNode["file_root"].as< string >( "vtkOutput" ) );
        outputs.push_back( sp );
      }
      else if( outputType == "restart" )
      {
        auto sp = std::make_shared< RestartOutput >( i,
                                                     subNode["every"].as< string >( "" ),
                                                     subNode["at"].as< std::vector< string > >( std::vector< string >() ) );
        outputs.push_back( sp );
      }
      else
      {
        GEOS_WARNING( "Discarded output \"" << outputType << "\"" );
      }
    }
  }
}

void operator>>( const YAML::Node & node,
                 Deck & deck )
{
  Simulation simulation;
  node["simulation"] >> simulation;
  deck.setSimulation( simulation );

  std::vector< std::shared_ptr< Output > > outputs;
  node["outputs"] >> outputs;
  deck.setOutputs( outputs );
}

void convert( string const & stableInputFileName, xmlWrapper::xmlDocument & doc )
{
  xml_node problem = doc.appendChild( "Problem" );

  YAML::Node const input = YAML::LoadFile( stableInputFileName );
  Deck deck;
  input >> deck;

  deck.fillProblemXmlNode( problem );

  doc.getPugiDocument().save( std::cout, "    ", pugi::format_indent | pugi::format_indent_attributes );
}

}
