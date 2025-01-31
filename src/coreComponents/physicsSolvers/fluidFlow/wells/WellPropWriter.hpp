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
 * @file WellPropWriter.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPROPWRITER_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPROPWRITER_HPP

#include <map>

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"

#include "mesh/PerforationFields.hpp"

namespace geos
{

/******************************** WellPropWriter ********************************/

class WellPropWriter
{
public:
  template< bool cond, typename U >
  using resolvedType  = typename std::enable_if< cond, U >::type;

  struct prop_wrapper
  {
    virtual ~prop_wrapper() = default;
  };

  template< typename T >
  struct typed_prop_wrapper : public prop_wrapper
  {
    typed_prop_wrapper( const T * prop ): m_prop( prop ){}
    const T * m_prop;
  };

  typedef std::multimap< std::type_index, prop_wrapper > type_info_to_prop_map;

  struct prop_writer
  {
    virtual void write_prop( const integer idx, std::ofstream & stream ) = 0;
  };

  struct perf_prop_writer
  {
    virtual void write_prop( const integer j, const integer er, const integer esr, const integer ei, std::ofstream & stream ) = 0;
  };
  template< typename T > struct typed_1d_prop_writer : public prop_writer
  {
    typedef T value_type;
    typedef typed_prop_wrapper< value_type > prop_wrapper_type;

    typed_1d_prop_writer( value_type & prop ): m_prop( prop ){}
    value_type m_prop;

    virtual void write_prop( const integer idx, std::ofstream & stream )
    {
      stream << "," <<  (m_prop)[idx];
    }
  };
  template< typename T > struct typed_1d_perf_res_prop_writer : public perf_prop_writer
  {
    typedef T value_type;
    typedef typed_prop_wrapper< value_type > prop_wrapper_type;

    typed_1d_perf_res_prop_writer( value_type const & prop ): m_prop( prop ){}
    value_type const m_prop;

    virtual void write_prop( const integer g, const integer er, const integer esr, const integer ei, std::ofstream & stream )
    {
      GEOS_UNUSED_VAR( g );
      stream << "," <<  m_prop[er][esr][ei];
    }
  };
  template< typename T > struct typed_2d_perf_res_prop_writer1 : public perf_prop_writer
  {
    typedef T value_type;
    typedef typed_prop_wrapper< value_type > prop_wrapper_type;

    typed_2d_perf_res_prop_writer1( value_type & prop, integer dim2 ): m_prop( prop ), m_Dim2( dim2 ){}
    value_type const m_prop;
    const integer m_Dim2;

    virtual void write_prop( const integer j, const integer er, const integer esr, const integer ei, std::ofstream & stream )
    {
      GEOS_UNUSED_VAR( er );
      GEOS_UNUSED_VAR( esr );
      GEOS_UNUSED_VAR( ei );
      for( integer i=0; i<m_Dim2; i++ )
        stream << "," <<  m_prop[j][i];
    }
  };

  template< typename T > struct typed_2d_perf_res_prop_writerf : public perf_prop_writer
  {
    typedef T value_type;
    typedef typed_prop_wrapper< value_type > prop_wrapper_type;

    typed_2d_perf_res_prop_writerf( value_type & prop, integer dim2 ): m_prop( prop ), m_Dim2( dim2 ){}
    value_type const m_prop;
    const integer m_Dim2;

    virtual void write_prop( const integer g, const integer er, const integer esr, const integer ei, std::ofstream & stream )
    {
      GEOS_UNUSED_VAR( g );
      for( integer i=0; i<m_Dim2; i++ )
        stream << "," <<  m_prop[er][esr][ei][0][i];
    }
  };
  template< typename T > struct typed_2d_perf_res_prop_writer : public perf_prop_writer
  {
    typedef T value_type;
    typedef typed_prop_wrapper< value_type > prop_wrapper_type;

    typed_2d_perf_res_prop_writer( value_type & prop, integer dim2 ): m_prop( prop ), m_Dim2( dim2 ){}
    value_type const m_prop;
    const integer m_Dim2;

    virtual void write_prop( const integer g, const integer er, const integer esr, const integer ei, std::ofstream & stream )
    {
      GEOS_UNUSED_VAR( g );
      for( integer i=0; i<m_Dim2; i++ )
        stream << "," <<  m_prop[er][esr][ei][i];
    }
  };
  template< typename T > struct typed_2d_prop_writer : public prop_writer
  {
    typedef T value_type;
    typedef typed_prop_wrapper< value_type > prop_wrapper_type;

    typed_2d_prop_writer( value_type & prop, integer dim2 ): m_prop( prop ), m_Dim2( dim2 ){}
    value_type m_prop;
    const integer m_Dim2;

    virtual void write_prop( const integer idx, std::ofstream & stream )
    {

      for( integer i=0; i<m_Dim2; i++ )
        stream << "," <<  m_prop[idx][i];
    }
  };
  template< typename T > struct typed_f2d_prop_writer : public prop_writer
  {
    typedef T value_type;
    typedef typed_prop_wrapper< value_type > prop_wrapper_type;

    typed_f2d_prop_writer( value_type & prop, integer dim2 ): m_prop( prop ), m_Dim2( dim2 ){}
    value_type m_prop;
    const integer m_Dim2;

    virtual void write_prop( const integer idx, std::ofstream & stream )
    {
      for( integer i=0; i<m_Dim2; i++ )
        stream << "," <<  m_prop[idx][0][i];
    }
  };

  template< typename T > struct typed_f2d_perf_prop_writer : public perf_prop_writer
  {
    typedef T value_type;
    typedef typed_prop_wrapper< value_type > prop_wrapper_type;

    typed_f2d_perf_prop_writer( value_type & prop, integer dim2 ): m_prop( prop ), m_Dim2( dim2 ){}
    value_type m_prop;
    const integer m_Dim2;

    virtual void write_prop( const integer g, const integer er, const integer esr, const integer ei, std::ofstream & stream )
    {
      GEOS_UNUSED_VAR( g );
      for( integer i=0; i<m_Dim2; i++ )
        stream << "," <<  m_prop[er][esr][ei][0][i];
    }
  };
  template< typename T > struct typed_3d_prop_writer : public prop_writer
  {
    typedef T value_type;
    typedef typed_prop_wrapper< value_type > prop_wrapper_type;

    typed_3d_prop_writer( value_type & prop, integer dim2, integer dim3 ): m_prop( prop ), m_Dim2( dim2 ), m_Dim3( dim3 ){}
    value_type m_prop;
    const integer m_Dim2;
    const integer m_Dim3;

    virtual void write_prop( const integer idx, std::ofstream & stream )
    {
      for( integer i=0; i<m_Dim2; i++ )
        for( integer j=0; j<m_Dim3; j++ )
          stream << "," <<  m_prop[idx][i][j];
    }
  };

  template< typename T > struct typed_f3d_prop_writer : public prop_writer
  {
    typedef T value_type;
    typedef typed_prop_wrapper< value_type > prop_wrapper_type;

    typed_f3d_prop_writer( value_type & prop, integer dim2, integer dim3 ): m_prop( prop ), m_Dim2( dim2 ), m_Dim3( dim3 ){}
    value_type m_prop;
    const integer m_Dim2;
    const integer m_Dim3;

    virtual void write_prop( const integer idx, std::ofstream & stream )
    {
      for( integer i=0; i<m_Dim2; i++ )
        for( integer j=0; j<m_Dim3; j++ )
          stream << "," <<  m_prop[idx][0][i][j];
    }
  };

  template< typename T > struct typed_f3d_perf_res_prop_writer : public perf_prop_writer
  {
    typedef T value_type;
    typedef typed_prop_wrapper< value_type > prop_wrapper_type;

    typed_f3d_perf_res_prop_writer( value_type & prop, integer dim2, integer dim3 ): m_prop( prop ), m_Dim2( dim2 ), m_Dim3( dim3 ){}
    value_type m_prop;
    const integer m_Dim2;
    const integer m_Dim3;

    virtual void write_prop( const integer g, const integer er, const integer esr, const integer ei, std::ofstream & stream )
    {
      GEOS_UNUSED_VAR( g );
      for( integer i=0; i<m_Dim2; i++ )
        for( integer j=0; j<m_Dim3; j++ )
          stream << "," <<  m_prop[er][esr][ei][0][i][j];
    }
  };

  typedef std::unordered_map< std::string, prop_wrapper > prop_to_writer_map;
  typedef std::vector< prop_writer * > prop_to_writer_vec;
  typedef std::vector< perf_prop_writer * > perfprop_to_writer_vec;

  WellPropWriter ()
  {
    m_initialized=0;
  }
  WellPropWriter( WellPropWriter && ) = default;
  ~WellPropWriter()
  {
    if( m_outputFile.is_open() )
    {
      m_outputFile.close();
    }
    if( m_perfOutputFile.is_open() )
    {
      m_perfOutputFile.close();
    }
  }
  void initialize_perf( int myrank, const string & outputDir, const string & wellName, PerforationData const & perfData )
  {
    m_perfOutputFile.open( outputDir + "/" + wellName + "_perf_"+ std::to_string( myrank )+".csv" );
    m_numPerforations=perfData.size();
    m_resElementRegion = perfData.getField< fields::perforation::reservoirElementRegion >();
    m_resElementSubRegion = perfData.getField< fields::perforation::reservoirElementSubRegion >();
    m_resElementIndex = perfData.getField< fields::perforation::reservoirElementIndex >();
    m_perfWellElemIndex = perfData.getField< fields::perforation::wellElementIndex >();
    m_perfResElemGlobalIndex = perfData.getField< fields::perforation::reservoirElementGlobalIndex >();
  }
  void initialize_seg( int myrank, const string & outputDir, const string & wellName,
                       const arrayView1d< string const > & phaseNames,
                       const arrayView1d< string const > & componentNames,
                       WellElementSubRegion & subRegion )
  {
    m_outputFile.open( outputDir + "/" + wellName + "_seg_" + std::to_string( myrank )+".csv" );
    m_numSegments = subRegion.size();
    m_phaseNames = phaseNames;
    m_componentNames = componentNames;
    m_numPhase = m_phaseNames.size();
    m_numComponent = m_componentNames.size();
    m_elemGhostRank = subRegion.ghostRank();
    m_globalWellElementIndex = subRegion.getGlobalWellElementIndex();
  }

  template< typename T >
  void registerSegProp( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      m_header.push_back( name );
    m_propWriterVec.push_back( new typed_1d_prop_writer( prop ));
  }
  template< typename T >
  void registerPerfResProp( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      m_perfHeader.push_back( name );
    m_perfPropWriterVec.push_back( new typed_1d_perf_res_prop_writer( prop ));
  }
  template< typename T >
  void registerPerfComponentProp( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & c : m_componentNames )
      {
        m_perfHeader.push_back( name+"_"+c );
      }
    m_perfPropWriterVec.push_back( new typed_2d_perf_res_prop_writer1( prop, m_numComponent ));
  }
  template< typename T >
  void registerSeg2dProp( std::vector< std::string > const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & c : name )
      {
        m_header.push_back( c );
      }
    m_propWriterVec.push_back( new typed_2d_prop_writer( prop, name.size() ));
  }
  template< typename T >
  void registerPerf2dProp( std::vector< std::string > const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & c : name )
      {
        m_perfHeader.push_back( c );
      }
    m_perfPropWriterVec.push_back( new typed_2d_perf_res_prop_writer1( prop, name.size() ));
  }
  template< typename T >
  void registerPerfResComponentProp( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & c : m_componentNames )
      {
        m_perfHeader.push_back( name+"_"+c );
      }
    m_perfPropWriterVec.push_back( new typed_2d_perf_res_prop_writer( prop, m_numComponent ));
  }
  template< typename T >
  void registerPerfResPhaseComponentProp( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & p : m_phaseNames )
      {
        for( const auto & c : m_componentNames )
        {
          m_perfHeader.push_back( name+"_"+p+"_"+c );
        }
      }
    m_perfPropWriterVec.push_back( new typed_f3d_perf_res_prop_writer( prop, m_numPhase, m_numComponent ));
  }
  template< typename T >
  void registerSegPhasePropf( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & p : m_phaseNames )
      {
        m_header.push_back( name+"_"+p );
      }
    m_propWriterVec.push_back( new typed_f2d_prop_writer( prop, m_numPhase ));
  }
  template< typename T >
  void registerSegPhasePropDerf( std::string const & name, const std::vector< string > & der_names, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & p : m_phaseNames )
      {
        for( const auto & d:der_names )
          m_header.push_back( name+"_"+p +"_"+d );
      }

    m_propWriterVec.push_back( new typed_f3d_prop_writer( prop, m_numPhase, der_names.size()));

  }
  template< typename T >
  void registerSegPhasePropDer( std::string const & name, const std::vector< string > & der_names, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & p : m_phaseNames )
      {
        for( const auto & d:der_names )
          m_header.push_back( name+"_"+p +"_"+d );
      }
    m_propWriterVec.push_back( new typed_3d_prop_writer( prop, m_numPhase, der_names.size()));
  }
  template< typename T >
  void registerSegPhaseProp( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & p : m_phaseNames )
      {
        m_header.push_back( name+"_"+p );
      }

    m_propWriterVec.push_back( new typed_2d_prop_writer( prop, m_numPhase ));
  }
  template< typename T >
  void registerPerfResPhaseProp( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & p : m_phaseNames )
      {
        m_perfHeader.push_back( name+"_"+p );
      }

    m_perfPropWriterVec.push_back( new typed_2d_perf_res_prop_writer( prop, m_numPhase ));
  }
  template< typename T >
  void registerPerfResPhasePropf( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & p : m_phaseNames )
      {
        m_perfHeader.push_back( name+"_"+p );
      }

    m_perfPropWriterVec.push_back( new typed_2d_perf_res_prop_writerf( prop, m_numPhase ));
  }
  template< typename T >
  void registerPerfPhasePropf( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & p : m_phaseNames )
      {
        m_perfHeader.push_back( name+"_"+p );
      }

    m_perfPropWriterVec.push_back( new typed_f2d_perf_prop_writer( prop, m_numPhase ));
  }
  template< typename T >
  void registerPerfComponentPropf( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & p : m_componentNames )
      {
        m_perfHeader.push_back( name+"_"+p );
      }

    m_perfPropWriterVec.push_back( new typed_f2d_perf_prop_writer( prop, m_numComponent ));
  }

  template< typename T >
  void registerSegComponentProp( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & c : m_componentNames )
      {
        m_header.push_back( name+"_"+c );
      }
    m_propWriterVec.push_back( new typed_2d_prop_writer( prop, m_numComponent ));
  }
  template< typename T >
  void registerSegPhaseComponentProp( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & p : m_phaseNames )
      {
        for( const auto & c : m_componentNames )
        {
          m_header.push_back( name+"_"+p+"_"+c );
        }
      }
    m_propWriterVec.push_back( new typed_3d_prop_writer( prop, m_numPhase, m_numComponent ));
  }
  template< typename T >
  void registerSegPhaseComponentPropf( std::string const & name, const T & prop )
  {
    if( m_initialized == 0 )
      for( const auto & p : m_phaseNames )
      {
        for( const auto & c : m_componentNames )
        {
          m_header.push_back( name+"_"+p+"_"+c );
        }
      }
    m_propWriterVec.push_back( new typed_f3d_prop_writer( prop, m_numPhase, m_numComponent ));
  }
  void write( real64 time, real64 dt, integer cycle, integer subevent, integer timeStep, integer newtonIter, integer numTimeStepCuts )
  {
    if( m_initialized == 0 )
    {
      m_outputFile << "Time,Dt,Cycle,SubEvent,TimeStep,NewtonIteration,TimeStepCuts,Element";
      for( auto i :m_header )
      {
        m_outputFile << "," << i;
      }
      m_outputFile << std::endl;
      // for perf data
      m_perfOutputFile << "Time,Dt,Cycle,SubEvent,TimeStep,NewtonIteration,TimeStepCuts,ResElement,WellElement";
      for( auto i :m_perfHeader )
      {
        m_perfOutputFile << "," << i;
      }
      m_perfOutputFile << std::endl;
      m_initialized = 1;
    }

    for( integer j =0; j < m_numSegments; j++ )
    {
      if( m_elemGhostRank[j] < 0 )
      {
        m_outputFile <<   time << "," <<  dt << "," <<  cycle << "," <<  subevent << "," << timeStep << "," << newtonIter << "," << numTimeStepCuts<<","<<m_globalWellElementIndex[j];
        for( auto i : m_propWriterVec )
        {
          i->write_prop( j, m_outputFile );
        }
        m_outputFile << std::endl;
      }
    }

    for( integer j =0; j < m_numPerforations; j++ )
    {
      if( m_elemGhostRank[j] < 0 )
      {
        localIndex const er  = m_resElementRegion[j];
        localIndex const esr = m_resElementSubRegion[j];
        localIndex const ei  = m_resElementIndex[j];
        localIndex const iwelem = m_perfWellElemIndex[j];
        m_perfOutputFile << time << "," <<  dt << "," <<  cycle << "," <<  subevent << "," << timeStep << "," << newtonIter << "," << numTimeStepCuts<<","<<m_perfResElemGlobalIndex[j] << "," <<
          m_globalWellElementIndex[iwelem];
        for( auto i : m_perfPropWriterVec )
        {
          i->write_prop( j, er, esr, ei, m_perfOutputFile );
        }
        m_perfOutputFile << std::endl;
      }
    }
    m_propWriterVec.clear();
    m_perfPropWriterVec.clear();
  }

protected:
  integer m_initialized;
  arrayView1d< string const > m_phaseNames;
  arrayView1d< string const > m_componentNames;
  integer m_numSegments;
  integer m_numComponent;
  integer m_numPhase;

  std::vector< string > m_header;
  std::ofstream m_outputFile;
  prop_to_writer_vec m_propWriterVec;

  // perforation data
  integer m_numPerforations;
  arrayView1d< localIndex const >  m_resElementRegion;
  arrayView1d< localIndex const >  m_resElementSubRegion;
  arrayView1d< localIndex const >  m_resElementIndex;
  arrayView1d< localIndex const >  m_perfWellElemIndex;
  arrayView1d< globalIndex const >  m_perfResElemGlobalIndex;

  std::vector< string > m_perfHeader;
  std::ofstream m_perfOutputFile;
  perfprop_to_writer_vec m_perfPropWriterVec;

  /// Global index of local element
  arrayView1d< globalIndex const >  m_globalWellElementIndex;
  arrayView1d< integer const >  m_elemGhostRank;
};


} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_PERFORATIONFLUXLKERNELS_HPP
