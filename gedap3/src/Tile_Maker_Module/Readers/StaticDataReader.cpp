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
 * StaticDataReader.cpp
 *
 *  Created on: 13 Oct 2015
 *      Author: alix
 */

#include "StaticDataReader.h"

using std::string;
using std::shared_ptr;
using std::make_shared;
using std::cout;
using std::endl;
using std::get;

using netCDF::NcFile;
using netCDF::exceptions::NcException;

StaticDataReader::StaticDataReader(const string latlon_filename, const string land_sea_mask_filename,
								   const string elevation_filename,
								   shared_ptr<AbstractLogger> logger)
    : AbstractReader(logger)
{
	// Try to open input file
	isReady = true;

	try {
        logger->write_log(LogLevel::_DEBUG_, "StaticDataReader", "StaticDataReader",
                          "Opening NetCDF LatLonFile file " + latlon_filename);
		LatLonFile = make_shared<NcFile>(latlon_filename, NcFile::read);     

	}
	catch (NcException &c) {
		/* NetCDF-CXX 4.3 required
		logger->write_log(LogLevel::_ERROR_, "StaticDataReader", "StaticDataReader",
                          "Error while attempting to open lat/lon file " + latlon_filename +
                          " [NetCDF err code " + to_string(c.errorCode()) + "]");*/
		logger->write_log(LogLevel::_ERROR_, "StaticDataReader", "StaticDataReader",
						  "Error while attempting to open lat/lon file " + latlon_filename);
		isReady = false;
	}

	try {
        logger->write_log(LogLevel::_DEBUG_, "StaticDataReader", "StaticDataReader",
                          "Opening NetCDF LandSeaMaskFile file " + land_sea_mask_filename);
		LandSeaMaskFile = make_shared<NcFile>(land_sea_mask_filename, NcFile::read);

	}
	catch (NcException &c) {
		/* NetCDF-CXX 4.3 required
		logger->write_log(LogLevel::_ERROR_, "StaticDataReader", "StaticDataReader",
                          "Error while attempting to open land sea mask file " + land_sea_mask_filename +
                          " [NetCDF err code " + to_string(c.errorCode()) + "]");*/
		logger->write_log(LogLevel::_ERROR_, "StaticDataReader", "StaticDataReader",
						  "Error while attempting to open land sea mask file " + land_sea_mask_filename);
        isReady = false;
    }

    try {
        ElevationFile = make_shared<NcFile>(elevation_filename, NcFile::read);
    }
    catch (NcException &c) {
		/* NetCDF-CXX 4.3 required
        logger->write_log(LogLevel::_ERROR_, "StaticDataReader", "StaticDataReader",
                          "Error while attempting to open elevation file " + elevation_filename +
                          " [NetCDF err code " + to_string(c.errorCode()) + "]");*/
		logger->write_log(LogLevel::_ERROR_, "StaticDataReader", "StaticDataReader",
						  "Error while attempting to open elevation file " + elevation_filename);
		isReady = false;
	}


    if (!isReady) {
        logger->write_log(LogLevel::_WARNING_, "StaticDataReader", "StaticDataReader",
                          "One of the static data files could not be opened");
    }
    else {
        logger->write_log(LogLevel::_INFO_, "StaticDataReader", "StaticDataqReader",
                          "Creation of StaticDataReader instance successful");
    }
}


void StaticDataReader::extract_data(StaticData& static_data, TileDimensions& tile_dims)
{
    tile_dims.set_start_count_vectors(start, count);
    read_latitude_and_longitude(static_data);
	read_land_sea_mask(static_data);
	read_elevation(static_data);
}


/* ****************************************************************************************************	*
 * 											PRIVATE METHODS												*
 * ****************************************************************************************************	*/

void StaticDataReader::read_latitude_and_longitude(StaticData& static_data)
{
	try{
		LatLonFile->getVar("latitude").getVar(start, count, static_data.latitude.data());
		LatLonFile->getVar("longitude").getVar(start, count, static_data.longitude.data());
	}
	catch (NcException& c){
		cout << "error when reading lat/lon file" << endl;
   }
}


void StaticDataReader::read_land_sea_mask(StaticData& static_data)
{
	try{
		LandSeaMaskFile->getVar("land_sea_mask").getVar(start, count, static_data.landCover.data());
	}

	catch (NcException& c){
		cout << "error when reading land sea mask file" << endl;
   }
}


void StaticDataReader::read_elevation(StaticData& static_data)
{
	try{
		ElevationFile->getVar("elevation").getVar(start, count, static_data.elevation.data());
	}
	catch (NcException& c){
		cout << "error when reading elevation file" << endl;
   }
}



