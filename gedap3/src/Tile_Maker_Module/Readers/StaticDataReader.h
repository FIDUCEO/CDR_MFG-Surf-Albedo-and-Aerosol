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
 * StaticDataReader.h
 *
 *  Created on: 13 Oct 2015
 *      Author: alix
 */

#ifndef TILE_MAKER_MODULE_STATICDATAREADER_H_
#define TILE_MAKER_MODULE_STATICDATAREADER_H_

#include <netcdf>

#include "AbstractReader.h"
#include "Satellite_Parameters.h"

//#define DEBUG_LAT_LON
//#define DEBUG_LAND_SEA_MASK
//#define DEBUG_ELEVATION
//#define DEBUG_AEROSOLS_ELEVATION

/**
 * \ingroup Tile_Maker
 *
 * \brief   Read longitude, latitude, land cover, elevation and aerosol layer height from external files and
 *          copy the data in an instance of StaticData structure for each geographical tile.
 * \details \sa StaticData
 *
 * \warning External data are assumed to be already in the right projection.\n
 *          The static data (latitude, longitude, land cover, elevation) are projected on the SEVIRI grid using the program located in:\n
 *          svn://rayser01/RAYFERENCE/Algorithms/GEDAP/SEVIRI_Land_Cover_and_Elevation
 */
class StaticDataReader : public AbstractReader
{
	// input netCDF files
	std::shared_ptr<netCDF::NcFile> LatLonFile;
	std::shared_ptr<netCDF::NcFile> LandSeaMaskFile;
	std::shared_ptr<netCDF::NcFile> ElevationFile;

	//! To read a block of a dataset
	std::vector<size_t> start;
	//! To read a block of a dataset
	std::vector<size_t> count;

public:
	StaticDataReader(const std::string latlon_filename, const std::string land_sea_mask_filename,
					 const std::string elevation_filename, std::shared_ptr<AbstractLogger> logger);

    virtual ~StaticDataReader() = default;
//    virtual ~StaticDataReader() {  logger->write_log(LogLevel::_DEBUG_, "StaticDataReader", "~StaticDataReader",
//                                                           "Calling destructor StaticDataReader");}


	/**
	 * \brief For a given tile, extract longitude, latitude, land cover, elevation and aerosol layer height from external files
	 *        and copy it into an instance of StaticData structure.
	 * \details \sa StaticData
	 */
    void extract_data(StaticData& static_data, TileDimensions& tile_dims);

private:
	void read_latitude_and_longitude(StaticData& static_data);

	void read_land_sea_mask(StaticData& static_data);

	void read_elevation(StaticData& static_data);

};

#endif /* TILE_MAKER_MODULE_STATICDATAREADER_H_ */
