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
 * MVIRITileMaker.cpp
 *
 *  Created on: May 2017
 *      Author: Cedric
 */

#include "MFGTileMaker.h"
#include "../Tile_Maker_Module/Readers/MFGObservationsReader.h"
#include "../Tile_Maker_Module/Readers/MFGMaskReader.h"

using std::ifstream;
using std::string;
using std::vector;
using std::map;
using std::shared_ptr;
using std::make_shared;

using netCDF::NcFile;
using netCDF::exceptions::NcException;

MFGTileMaker::MFGTileMaker(const vector<string> &v_bands, const vector<TileDimensions_Tuple> &v_tiles_as_tuples,
                               shared_ptr<AbstractLogger> logger)
        : TileMaker(v_bands, v_tiles_as_tuples, logger) {
}

void MFGTileMaker::make_input_tile(map<string, string> m_input_data_files, const string date, const string time,
                                   const int dayOfMonth, const int counter, const int slot,
                                   const string path_image_file, const string path_mask_file,
                                   const string path_modparams_file, const string path_elevation_file,
                                   const string path_geopotential_file, const string aerosol_height_file_1,
                                   const string aerosol_height_file_2, const int reliability_cfc_retrieval, const bool heightIsConstant,
                                   const bool useCMASKSCORE, const bool usePCFC,
                                   const bool ecmwf_constant, const bool cold_start, const bool daily_cloudmaks) {

    clock_t start;
    double duration_read_observation = 0.0;
    double duration_read_era_interim = 0.0;
    double duration_read_mask = 0.0;
    double duration_extract_observation = 0.0;
    double duration_extract_era_interim = 0.0;
    double duration_extract_mask = 0.0;
    double duration_write_tile = 0.0;

    string input_tile_file = "";
    bool b_tile_to_create = false;

#ifdef DEBUG_TILE_MAKER
    cout << "map m_input_data_files: " << endl;
    for (pair<string, string> item : m_input_data_files)
        cout << item.first << " " << item.second << endl;
    cout << endl;
#endif

    // Check if there is Tiles to create
    for (TileDimensions &tile : v_tiles) {
        input_tile_file = m_input_data_files.at(tile.st_ID);
        if (!ifstream((input_tile_file).c_str()).good()) {
            b_tile_to_create = true;
        }
        else if (ifstream((input_tile_file).c_str()).good() && cold_start) {
            remove((input_tile_file).c_str());
            b_tile_to_create = true;
        }
    }
    // Create an instance of MFGObservationsReader, ModelParamsReader and MaskReader
    // TODO: check if shared_ptr is appropriate (maybe unique_ptr or just instantiating the object would be better)

    setup(path_mask_file, logger);


    if (b_tile_to_create) {
        logger->write_log(LogLevel::_INFO_, "MFGTileMaker", "make_input_tile",
                          "Making tiles associated with date-time " + date + " : " + time);

        mp_observations_reader->setup(path_image_file, logger);
        mp_mask_reader->setup(path_mask_file, reliability_cfc_retrieval, logger);

        start = clock();

        if (!ecmwf_constant) {
            int ecmwf_image_slot = 0.0;
            if (slot < 12) {
             ecmwf_image_slot = (dayOfMonth-1)*4 ;
            } else if (slot < 24) {
              ecmwf_image_slot = (dayOfMonth-1)*4 + 1;
            } else if (slot < 36) {
              ecmwf_image_slot = (dayOfMonth-1)*4 + 2 ;
            } else if (slot < 48) {
                ecmwf_image_slot = (dayOfMonth-1)*4 + 3 ;
              }

// The switch between nextData and currentData should not happen at every loop
// but only at the beginning of the tile creation (first tile, first processed slot) and at each multiple of 12
// to make sure this step was performed at the very beginning, an additional argument has been passed to the tile_maker, i.e. the index of the period

            if (counter == 0 || slot % 12 == 0){
               m_modelParamsReader->compute_gradient_for_interpolation(path_modparams_file, ecmwf_image_slot+1);
            }

            m_modelParamsReader->create_full_era_interim_image_at_step(slot);
            duration_read_era_interim += (double) (clock() - start);
        }

        // read "full image" mask data

        start = clock();

        if (daily_cloudmaks){
            mp_mask_reader->IODC_read(slot);
        }
        else {
            mp_mask_reader->read();
        }

        duration_read_mask += (double) (clock() - start);

        start = clock();
        for (string st_channel : v_bands)
            mp_observations_reader->add_channel(st_channel);

        if (mp_observations_reader->is_ready() && mp_mask_reader->is_ready()) {
            /* ****************************************	*
             *          SPATIAL LOOP (TILES)          *
             * ****************************************	*/
            mp_observations_reader->read();
            duration_read_observation += (double) (clock() - start);

            for (TileDimensions &tile : v_tiles) {

                // Get tile filename
                input_tile_file = m_input_data_files.at(tile.st_ID);
//                static_tile_file = m_static_data_files.at(tile.st_ID);

                if (!ifstream((input_tile_file).c_str()).good()) {
                    // Create an instance of InputData
                    InputData input_data(tile.nb_lines, tile.nb_columns, v_bands);

                    // Extract Observations
                    start = clock();
                    mp_observations_reader->extract_data(input_data, tile);

                    duration_extract_observation += (double) (clock() - start);

                    // Extract Model Parameters
                    try {

                        start = clock();
                        m_modelParamsReader->extract_data(input_data, path_elevation_file, path_geopotential_file,
                                                          aerosol_height_file_1, aerosol_height_file_2, dayOfMonth,
                                                          heightIsConstant, tile, slot % 12, ecmwf_constant);

                        duration_extract_era_interim += (double) (clock() - start);
                    }
                    catch (NcException &c) {
                        string error_message =
                                "Error while extracting model parameters: aborting input tile generation";
                        logger->write_log(LogLevel::_CRITICAL_, "TileMaker", "make_input_tile", error_message);
                        throw std::runtime_error(error_message);
                    }

                    // Extract Mask
                    start = clock();

                    mp_observations_reader->extract_data_quality_bitmask(input_data, tile);

                    mp_mask_reader->extract_data(input_data, useCMASKSCORE, usePCFC, tile);
                    duration_extract_mask += (double) (clock() - start);
                    // Save tile for the current tile
                    start = clock();

                    m_tile_netCDF_manager->save_input_data(input_tile_file, input_data);
                    duration_write_tile += (double) (clock() - start);

                } else {
                    logger->write_log(LogLevel::_INFO_, "TileMaker", "make_input_tile",
                                      "Input tile " + input_tile_file + " already exists.");
                }
            }
            if (!ecmwf_constant) {
                logger->write_log(LogLevel::_DEBUG_, "TileMaker", "make_input_tile",
                                  "time to read era-interim data: " +
                                  std::to_string(duration_read_era_interim / (double) CLOCKS_PER_SEC));
            }
            logger->write_log(LogLevel::_DEBUG_, "TileMaker", "make_input_tile",
                              "time to read observations: " + std::to_string(duration_read_observation / (double) CLOCKS_PER_SEC));
            logger->write_log(LogLevel::_DEBUG_, "TileMaker", "make_input_tile",
                              "time to read mask: " + std::to_string(duration_read_mask / (double) CLOCKS_PER_SEC));
            logger->write_log(LogLevel::_DEBUG_, "TileMaker", "make_input_tile",
                              "total time to read all data: " + std::to_string(
                                      (duration_read_observation + duration_read_era_interim + duration_read_mask) /
                                      (double) CLOCKS_PER_SEC));
            logger->write_log(LogLevel::_DEBUG_, "TileMaker", "make_input_tile",
                              "time to extract observations: " +
                              std::to_string(duration_extract_observation / (double) CLOCKS_PER_SEC));
            logger->write_log(LogLevel::_DEBUG_, "TileMaker", "make_input_tile",
                              "time to extract era-interim: " +
                              std::to_string(duration_extract_era_interim / (double) CLOCKS_PER_SEC));
            logger->write_log(LogLevel::_DEBUG_, "TileMaker", "make_input_tile",
                              "time to extract mask: " +
                              std::to_string(duration_extract_mask / (double) CLOCKS_PER_SEC));
            logger->write_log(LogLevel::_DEBUG_, "TileMaker", "make_input_tile",
                              "time to write tiles: " + std::to_string(duration_write_tile / (double) CLOCKS_PER_SEC));
            logger->write_log(LogLevel::_DEBUG_, "TileMaker", "make_input_tile",
                              "total time to extract all data: " + std::to_string(
                                      (duration_extract_observation + duration_extract_era_interim +
                                       duration_extract_mask +
                                       duration_write_tile) / (double) CLOCKS_PER_SEC));

            logger->write_log(LogLevel::_INFO_, "TileMaker", "make_input_tile",
                              "Generation of Input tiles associated with date-time " + date + " : " + time + " finished.");
        }
    }
    else { // if (b_tile_to_create)
        logger->write_log(LogLevel::_INFO_, "TileMaker", "make_input_tile",
                          "All Input tiles for date-time " + date + " : " + time +
                          " already exist.  There is nothing to do.");
    } // if (b_tile_to_create)
}

void MFGTileMaker::read_first_ERA_INTERIM(const string path_modparams_file, const int era_interim_slot) {
    m_modelParamsReader->read_first_model_params_data(path_modparams_file, era_interim_slot);
}


/* ****************************************************************************************************	*
 * 											PRIVATE METHODS												*
 * ****************************************************************************************************	*/

void MFGTileMaker::setup(const string inputFilename, shared_ptr<AbstractLogger> logger){
    //opening Angle image (Static EasyFCDR) here, so that it doesn't get destroyed every loop

    mp_observations_reader = make_shared<MFGObservationsReader>(logger);
    mp_mask_reader = make_shared<MFGMaskReader>(logger);

}
