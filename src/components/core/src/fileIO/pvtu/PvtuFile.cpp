/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file PvtuFile.cpp
 */

#include "PvtuFile.hpp"
#include "common/Logger.hpp"
#include <iostream>
#include <string.h>

#if USE_MPI
#include <mpi.h>
#endif

namespace geosx{

// PUBLIC METHODS
void PvtuFile::load( std::string const &filename) {
    int size = 0;
    int rank = 0;
#if USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    int nb_partition_per_core = 0;
    //int nb_total_partitions = 0;
    int remainder = 0;
    std::vector < std::string > children_files;



    // XML file parsing using pugixml tool
    pugi::xml_document pvtu_doc;
    pvtu_doc.load_file(filename.c_str());

    if( rank == 0) {
        check_xml_file_consistency(pvtu_doc);
    }
    vtu_files_list(pvtu_doc,children_files);
    // Retrieve the number of partitions
    int const nb_total_partitions = children_files.size();

    // Next part of this method is dedicated to the optimization of file loading
    //
    // IF nb_partitions > nb_mpi_process : Each process will load 
    // nb_partition / nb_mpi_process. Processes with a rank < nb_partition % nb_mpi_process
    // will load an additional partition
    //
    // IF nb_partitions < nb_mpi_process : the first nb_partitions process will
    // load ONE partition.
    if(nb_total_partitions > size) {
        if ( rank == 0 ) {
            std::cout << "WARNING : the number of partitions ("
                << nb_total_partitions <<") which will be loaded " 
                << "is greater that the number of processes on which GEOSX is running ("
                << size
                << "). Some processes will hold more than one partition !" << std::endl;
        }
        nb_partition_per_core = nb_total_partitions / size;
        remainder = nb_total_partitions % size;
    } else {
        nb_partition_per_core = 1;
        remainder = 0;
    }
    if( rank < remainder ) {
        vtu_file_names_.resize(nb_partition_per_core+1);
        vtu_files_.resize(nb_partition_per_core+1);
        for(int p_index = 0; p_index < nb_partition_per_core +1; ++p_index) {
            vtu_file_names_[p_index] =
                children_files[rank*(nb_partition_per_core+1) + p_index];
        }
    } else if( rank < nb_total_partitions) {
        vtu_file_names_.resize(nb_partition_per_core);
        vtu_files_.resize(nb_partition_per_core);
        for(int p_index = 0; p_index < nb_partition_per_core; ++p_index) {
            vtu_file_names_[p_index] =
                children_files[rank*(nb_partition_per_core) + p_index+remainder];
        }
    }
}

void PvtuFile::save( std::string const &fileName) {
    geos_abort("pvtu file save is not implemented yet");
}

void VtuFile::load( std::string const &fileName) {
}

void VtuFile::save( std::string const &fileName) {
    geos_abort("vtu file save is not implemented yet");
}

void VtuFile::check_xml_file_consistency(pugi::xml_document const & pvtu_doc) const{
}


//PRIVATE METHODS
void PvtuFile::check_xml_file_consistency(pugi::xml_document const & pvtu_doc) const {

    // VTKFile is the main node of the pvtufile
    auto const vtk_file =pvtu_doc.child("VTKFile");
    if( vtk_file.empty() ) {
        geos_abort("Main node VTKFile not found or empty in the parent file");
    }

    auto const pugrid = vtk_file.child("PUnstructuredGrid");
    if( pugrid.empty() ) {
        geos_abort("Node PUnstructuredGrid not found or empty in the parent file");
    }

    std::string const pugrid_child_names[4] = {"PPointData","PCellData","PPoints","Piece"};
    for( auto pugrid_child_name : pugrid_child_names ) {
        auto const pugrid_child = pugrid.child(pugrid_child_name.c_str());
        if(pugrid_child.empty()) {
            std::string const message = "Node " + pugrid_child_name +
                " not found or empty in the parent file";
            geos_abort(message);
        };
    }

    bool points_have_original_index_property = false;
    std::string const mandatory_attributes[4] =
    {"Name","type","format","NumberOfComponents"};
    for( auto ppdata_property : pugrid.child("PPointData").children() ) {
        if( ppdata_property.name() == static_cast< std::string >("DataArray")) {
            for( auto mandatory_attribute : mandatory_attributes ){
                auto const attribute =
                    ppdata_property.attribute( mandatory_attribute.c_str());
                if( attribute.empty() ) {
                    std::string const message ="Mandatory attribute " + mandatory_attribute +
                        " does not exist in a DataArray of PPointData";
                    geos_abort(message);
                }
            }
            if( ppdata_property.attribute("Name").value() == str_original_index_) {
                points_have_original_index_property = true;
            }
        }
    }
    if (!points_have_original_index_property) {
        geos_abort("Can't find any DataArray which contains the property "+
                str_original_index_ + " in PPointData");
    }

    bool cells_have_original_index_property = false;
    bool cells_have_partition_property = false;
    bool cells_have_region_property = false;
    for( auto ppdata_property : pugrid.child("PCellData").children() ) {
        if( ppdata_property.name() == static_cast< std::string >("DataArray")) {
            for( auto mandatory_attribute : mandatory_attributes ){
                auto const attribute =
                    ppdata_property.attribute( mandatory_attribute.c_str());
                if( attribute.empty() ) {
                    std::string const message ="Mandatory attribute " + mandatory_attribute +
                        " does not exist in a DataArray of PCellData";
                    geos_abort(message);
                }
            }
            if( ppdata_property.attribute("Name").value() == str_original_index_) {
                cells_have_original_index_property = true;
            }
            if( ppdata_property.attribute("Name").value() == str_partition_) {
                cells_have_partition_property = true;
            }
            if( ppdata_property.attribute("Name").value() == str_region_) {
                cells_have_region_property = true;
            }
        }
    }
    if (!cells_have_original_index_property) {
        geos_abort("Can't find any DataArray which contains the property "+
                str_original_index_ + " in PCellData");
    }
    if (!cells_have_region_property) {
        geos_abort("Can't find any DataArray which contains the property "+
                str_region_ + " in PCellData");
    }
    if (!cells_have_partition_property) {
        geos_abort("Can't find any DataArray which contains the property "+
                str_partition_ + " in PCellData");
    }

    bool ppoint_has_a_pdata_array = false;
    bool ppoint_has_a_pdata_array_with_points = false;
    for(auto ppoint_child : pugrid.child("PPoints").children() ) {
        if( ppoint_child.name() == static_cast< std::string >("PDataArray") ) {
            ppoint_has_a_pdata_array = true;
            if( ppoint_child.attribute("Name").as_string() ==
                    static_cast< std::string >("Points")) {
                ppoint_has_a_pdata_array_with_points = true;
                if (ppoint_child.attribute("NumberOfComponents").as_uint() != 3 ) {
                    geos_abort("GEOSX supports only 3D meshes");
                }
                break;
            }
        }
    }
    if( !ppoint_has_a_pdata_array ) {
        geos_abort("No PDataArray found in PPoints.");
    }
    if( !ppoint_has_a_pdata_array_with_points ) {
        geos_abort("No PDataArray named \"Points\" found");
    }

    for( auto pugrid_child : pugrid.children() ) {
        if(pugrid_child.name() == static_cast< std::string >("Piece")) {
            if( pugrid_child.attribute("Source").empty() ) {
                geos_abort("Piece nodes has to have an attribute Source not empty.");
            }
        }
    }
}
    void PvtuFile::vtu_files_list(
            pugi::xml_document const & pvtu_doc,
            std::vector < std::string > & vtu_files ) const{
        int nb_partitions = 0;
        for(auto child : pvtu_doc.child("VTKFile").child("PUnstructuredGrid").children()) {
            if( child.name() == static_cast< std::string > ("Piece") )
            {
                vtu_files.emplace_back(child.attribute("Source").as_string());
            }
        }
    }
}
