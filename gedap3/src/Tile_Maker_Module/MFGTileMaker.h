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
 * TileMaker.h
 *
 *  Created on: May 2017
 *      Author: Cedric
 */

#ifndef TILE_MAKER_MODULE_MFGTILEMAKER_H_
#define TILE_MAKER_MODULE_MFGTILEMAKER_H_

#include "TileMaker.h"
#include "Readers/MFGObservationsReader.h"
#include "Readers/MFGMaskReader.h"

//#define DEBUG_TILE_MAKER

/**
 * \defgroup Tile_Maker
 * \brief    Module responsible for the pre-processing, i.e. to aggregate external data and create the Static and Input netCDF files.
 * \warning  - The Reader classes assume that all the external data to be read are already in the MVIRI projection
 *           - The static data (latitude, longitude, land cover, elevation) are projected on the MVIRI grid using the program located in:\n
 *             svn://rayser01/RAYFERENCE/Algorithms/GEDAP/SEVIRI_Land_Cover_and_Elevation
 *           - The ECMWF data (TCO3, TCWV, surface pressure, wind speed, wind direction) are projected on the MVIRI grid using the program SOSIE located in:\n
 *             svn://rayser01/RAYFERENCE/Algorithms/GEDAP/SOSIE_ECMWF_interpolation
 *           - The external files containing the cloud mask are in the MVIRI projection.
 * \todo     Adapt all classes of the Tile_Maker module to be able to aggregate data for other satellites/radiometers.
 */


/**
 * \ingroup Tile_Maker
 *
 * \brief   Main class of the Tile_Maker module. It manages the creation of Static and Input Tile netCDF files.
 * \details This class is written as an interface dedicated to be used by the Python Scheduler module.
 */
class MFGTileMaker : public TileMaker {

private:
    std::shared_ptr<MFGObservationsReader> mp_observations_reader;
    std::shared_ptr<MFGMaskReader> mp_mask_reader;
    std::shared_ptr<ModelParamsReader> mp_modelParamsReader;

    int slot = 0.0;

    //std::string m_path_file_angles;

    void setup(const std::string inputFilename, std::shared_ptr<AbstractLogger> logger);

public:
    MFGTileMaker(const std::vector<std::string> &v_bands, const std::vector<TileDimensions_Tuple> &v_tiles_as_tuples,
                 std::shared_ptr<AbstractLogger> logger);

    virtual ~MFGTileMaker() = default;
//    virtual ~MFGTileMaker() {  logger->write_log(LogLevel::_DEBUG_, "MFGTileMaker", "~MFGTileMaker",
//                                                           "Calling destructor MFGTileMaker");}


    /*void make_static_data_tile(std::map<std::string, std::string> m_static_data_files, const std::string path_lat_lon_file,
                               const std::string path_land_sea_mask_file, const std::string path_elevation_file,
                               const bool cold_start) {
        TileMaker::make_static_data_tile(m_static_data_files, path_lat_lon_file, path_land_sea_mask_file,
                                         path_elevation_file, cold_start);
    }*/

    /**
     * \brief Generates the Input-Tile netCDF files for a given time slot.
     * \param m_input_data_files        map with items (key = tile ID, value = path of the Input-Tile netCDF file to create)
     * \param m_static_data_files       map with items (key = tile ID, value = path of the Static-Tile netCDF file)
     * \param date                      year, month and day written as YYYY-MM-DD (see Python code calling this method)
     * \param time                      hour and minutes written as HH-MM (see Python code calling this method)
     * \param dayOfYear                 day of the year [0 -> 365/366]
     * \param slot                      slot time of the day (required to interpolate the ECMWF data in time. MVIRI images = one image every 30 minutes = (24 * 2) slots)
     * \param path_image_file           path to the file containing the satellite image
     * \param path_mask_file            path to the file containing the cloud mask
     * \param path_next_modparams_file  path to the next ECMWF file
     * \param path_geopotential_field   path to the file containing the geopotential field
     */
    void make_input_tile(std::map<std::string, std::string> m_input_data_files, const std::string date, const std::string time,
                         const int dayOfMonth, const int counter, const int slot, const std::string path_image_file,
                         const std::string path_mask_file, const std::string path_modparams_file, const std::string path_elevation_file,
                         const std::string path_geopotential_file, const std::string aerosol_height_file_1,
                         const std::string aerosol_height_file_2, const int reliability_cfc_retrieval,
                         const bool heightIsConstant, const bool useCMASKSCORE, const bool usePCFC,
                         const bool ecmwf_constant, const bool cold_start, const bool daily_cloudmaks);

    void read_first_ERA_INTERIM(const std::string path_modparams_file, const int slot);

};

#endif /* TILE_MAKER_MODULE_MFGTILEMAKER_H_ */
