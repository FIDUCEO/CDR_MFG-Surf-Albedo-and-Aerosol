/**********************************************************************
* The copyrights for the GEDAP algorithm and computer codes, remain with 
* Rayference SPRL as an Intellectual Property Right 
*
* Any use, in source and binary forms, distribution, modification, for 
* commercial, business, research and any other purposes is
*
*                           *** STRICTLY FORBIDDEN ***
*
* GEDAP may not be distributed or sold to any other commercial, 
* business, research or other partners under any circumstances.
* In any case this comment is part of the code as Software Legal Information,
* it cannot be modified and must be kept as header for all source codes in 
* which it appears.
*
* All documents and software are protected under the international Copyright 
* Laws, and Rayference SPRL reserves all rights.
*
*
* Author       :     Rayference Copyright (c) 
*
*********************************************************************/
/*
 * TileMaker.cpp
 *
 *  Created on: 14 Oct 2015
 *      Author: alix
 */

#include <sstream>
#include <tuple>

#include "TileMaker.h"

using std::stringstream;
using std::ifstream;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::shared_ptr;
using std::make_shared;
using std::get;
using std::runtime_error;

TileMaker::TileMaker(const vector<string> &v_bands, const vector<TileDimensions_Tuple> &v_tiles_as_tuples,
					 shared_ptr<AbstractLogger> logger)
{
	this->logger = logger;

	this->v_bands = v_bands;
	v_tiles.reserve(v_tiles_as_tuples.size());

	for (const TileDimensions_Tuple &tuple_tile_dims : v_tiles_as_tuples) {
		v_tiles.push_back(TileDimensions(get<0>(tuple_tile_dims),
										 get<1>(tuple_tile_dims),
										 get<2>(tuple_tile_dims),
										 get<3>(tuple_tile_dims),
										 get<4>(tuple_tile_dims)));
	}

    stringstream ss;
    ss << "Channels: ";
    for (const string &channel : v_bands) ss << channel << " ";
    logger->write_log(LogLevel::_INFO_, "TileMaker", "TileMaker", ss.str());

    ss.clear();
    ss.str(string());
    ss << "Tile properties:" << endl;
    for (TileDimensions tile_dims : v_tiles) {
        ss << "    Index: " << tile_dims.st_ID << endl;
        ss << "        first line, last line, first column, last column: "
           << tile_dims.first_line << " "
           << tile_dims.last_line << " "
           << tile_dims.first_column << " "
           << tile_dims.last_column
           << endl;
        ss << "        nb lines, columns, pixels: "
           << tile_dims.nb_lines << " "
           << tile_dims.nb_columns << " "
           << tile_dims.nb_pixels;
    }
    logger->write_log(LogLevel::_INFO_, "TileMaker", "TileMaker", ss.str());

    m_tile_netCDF_manager = make_shared<TileNetCDFManager>(logger);
    m_modelParamsReader   = make_shared<ModelParamsReader>(logger);
}


void TileMaker::make_static_data_tile(map<string, string> m_static_data_files,
									  const string path_lat_lon_file, const string path_land_sea_mask_file,
									  const string path_elevation_file, const bool cold_start) {
	string static_tile_file = "";
	bool   b_tile_to_create = false;
	#ifdef DEBUG_TILE_MAKER
		cout << "map m_static_data_files: " << endl;
		for(pair<string, string> item : m_static_data_files)
			cout << item.first << " " << item.second << endl;
		cout << endl;
	#endif

	// Check if there is Tiles to create
	for (TileDimensions &tile : v_tiles) {
		static_tile_file = m_static_data_files.at(tile.st_ID);

		if (!ifstream((static_tile_file).c_str()).good()) {
            b_tile_to_create = true;
        }
		else if (ifstream((static_tile_file).c_str()).good() & cold_start) {
			remove((static_tile_file).c_str());
			b_tile_to_create = true;
		}
	}

    if (b_tile_to_create) {
        // Create an instance of StaticDataReader
        StaticDataReader staticDataReader(path_lat_lon_file, path_land_sea_mask_file, path_elevation_file, logger);

        if (!staticDataReader.is_ready()) {
            string error_message = "Error when looking up static data: aborting static tile generation";
            logger->write_log(LogLevel::_CRITICAL_, "TileMaker", "make_static_data_tile", error_message);
            throw runtime_error(error_message);
        } else {
            /* ****************************************	*
             *          SPATIAL LOOP (TILES)            *
             * ****************************************	*/
            for (TileDimensions &tile : v_tiles) {
                // Get static tile filename
                static_tile_file = m_static_data_files.at(tile.st_ID);
                if (!ifstream((static_tile_file).c_str()).good()) {
                    // Create an instance of Static Data
                    StaticData static_data(tile.nb_lines, tile.nb_columns);
                    // Extract data for the current tile + save the tile
                    staticDataReader.extract_data(static_data, tile);
                    // Save data into a netCDF file
                    m_tile_netCDF_manager->save_static_data(static_tile_file, static_data);
                } else {
                    logger->write_log(LogLevel::_INFO_, "TileMaker", "make_static_data_tile",
                                      "Static tile " + static_tile_file + " already exists.");
                }
            }

            logger->write_log(LogLevel::_INFO_, "TileMaker", "make_static_data_tile",
                              "Generation of Static tiles finished.");
        }
    } else {
        logger->write_log(LogLevel::_INFO_, "TileMaker", "make_static_data_tile",
                          "All Static tiles already exist. There is nothing to do.");
    }
}

//void TileMaker::read_first_ERA_INTERIM_file(const string path_modparams_file)
//{
//    // compute the first gradient for model parameters
//    m_modelParamsReader->read_first_model_params_data(path_modparams_file, slot);
//}
