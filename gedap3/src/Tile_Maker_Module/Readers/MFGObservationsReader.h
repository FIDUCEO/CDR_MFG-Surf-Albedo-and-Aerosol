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
 * MVIRIObservationsReader.h
 *
 *  Created on: 14 Oct 2015
 *      Author: alix
 */

#ifndef TILE_MAKER_MODULE_READERS_MFGOBSERVATIONSREADER_H_
#define TILE_MAKER_MODULE_READERS_MFGOBSERVATIONSREADER_H_

#include <iostream>
#include <cassert>
#include <netcdf>

#include "ObservationsReader.h"
#include "Satellite_Parameters.h"

//#define DEBUG_OBS
//#define DEBUG_COMMON_INPUT
//#define DEBUG_BRF
//#define DEBUG_COUNTS
//#define DEBUG_ANGLES
//#define DEBUG_ACQ_TIME
//#define DEBUG_RADIOMETRIC_NOISE

#define    NB_DETECTORS    42        // (3 * 11 + 9 = 42)

/**
 * \ingroup Tile_Maker
 *
 * \brief   Read counts, acquisition time and angles from external files and computes the TOA BRF
 *          (Top Of the Atmosphere Bid...) and radiometric noise.
 *          Extracted data are copied in an instance of Observations structure.
 * \details \sa Observations
 *
 * \warning External data are assumed to be already in the right projection.\n
 *
 * \note    We set a margin of 1 pixel around the tile to calculate derivatives as required to compute radiometric noise
 */
class MFGObservationsReader : public ObservationsReader {

private:
    // input files
    std::shared_ptr<netCDF::NcFile> mp_imageFile;

    //! index of the channel to be read. Starts at 0. HRVIS = 3
    short channel = -1;
    std::vector<std::string> v_channels;
    std::string st_channel_in_image_file;
    std::map<std::string, std::string> map_string_channel = {{"VIS", "01"}};

    //! We need 1 pixel margin to compute derivatives for radiometric noise
    const unsigned int margin = 1;

    //! Sun Zenith Angle [0 - 90]
    Eigen::ArrayXXf_RowMaj sza;
    //! Sun Azimuth Angle [0 - 360]
    Eigen::ArrayXXf_RowMaj saa;
    //! View Zenith Angle [0 - 90] (interpolated)
    Eigen::ArrayXXf_RowMaj vza;
    //! View Azimuth Angle [0 - 360] (interpolated)
    Eigen::ArrayXXf_RowMaj vaa;
    //! data quality bitmask
    Eigen::ArrayXXi_RowMaj dq_bitmask;


    //! Subsatellite point data
    struct {
        double lat_start;
        double lat_end;
        double lat_average;
        int lat_average_index; // Index of the average SSP in the latitude_values table
        double lon_start;
        double lon_end;
        double lon_average;
        int lon_average_index; // Index of the average SSP in the longitude_values table
    } ssp_info;
    //! List of latitude and longitude values in the angle table
    struct {
        Eigen::Array<double, Eigen::Dynamic, Eigen::RowMajor> latitude_values;
        Eigen::Array<double, Eigen::Dynamic, Eigen::RowMajor> longitude_values;
    } angle_table_info;

    //! Acquisition time
    std::vector<Eigen::ArrayXXf_RowMaj> acq_time;
    //! Additional info on acquisition time: min, max, average
    // WARNING: this data struct is incompatible with channel-specific data
    // (the channel dependency should be removed anyway from MVIRI-specific classes)
    struct {
        float min;
        float max;
        float mean;
    } acq_time_info;
    //! TOA BRF
    std::vector<Eigen::ArrayXXf_RowMaj> BRF;
    //! Post-processed uncertainty, aka total error covariance matrix (sigma_i^2 + sigma_r^2)
    std::vector<Eigen::ArrayXXf_RowMaj> radiometric_noise;


public:

    //! Construct from image and static files, and logger
    MFGObservationsReader(std::shared_ptr<AbstractLogger> logger);

    void setup(const std::string imageFilename, std::shared_ptr<AbstractLogger> logger);

     ~MFGObservationsReader() override = default;
//    ~MFGObservationsReader() override{ std::cout << "calling destructor observations reader" << std::endl;}


    //! Add an additional channel from std::string identifier
    void add_channel(const std::string st_channel) { v_channels.push_back(st_channel); }



    /**
     * \brief Read and store in buffers the data in the image file:
     * - acquisition time (and compute average time)
     * - angles (including satellite viewing angle table slice corresponding to the acquisition time frame)
     * - TOA BRF
     * - error data (also includes image-scale processing)
     */
    void read() override;

    /**
     * \brief For a given tile, extract and store the following data in an instance of the Observations structure:
     * - angles (and interpolation)
     * - acquisition time
     * - TOA BRF
     *
     *\details \sa Observations
     */
    void extract_data(Observations &obs_data, TileDimensions &tile_dims) override;

    //! Extract data_quality_bitmask
    void extract_data_quality_bitmask(Mask &mask_data, TileDimensions &tile_dims);


private:

    //! Read and process acquisition time info (esp. compute average acq. time)
    void read_acquisition_time();

    //! Read and process SSP info (esp. compute latitude and longitude of average SSP point)
    // Must be called after read_acquisition_time()
//    void read_ssp_info();

    //! Read angles
    // Must be called after read_ssp_info()
    void read_angles();

    //! Read TOA BRF
    void read_toa_brf();

    //! Read and compute uncertainty data
    // Must be called after read_toa_brf()
    void read_error();

    //! Extract angles (SZA, SAA, VZA and VAA); perform interpolation in the case of satellite angles
    void extract_angles(Observations &obs_data, TileDimensions &tile_dims);

    //! Extract channel wavelength, acquisition time, toa brf and radiometric noise data
    void extract_channel_observations(BandObservations &band_obs_data, TileDimensions &tile_dims);

    //! Read data_quality_bitmask
    void read_data_quality_bitmask();



};

#endif /* TILE_MAKER_MODULE_READERS_MFGOBSERVATIONSREADER_H_ */
