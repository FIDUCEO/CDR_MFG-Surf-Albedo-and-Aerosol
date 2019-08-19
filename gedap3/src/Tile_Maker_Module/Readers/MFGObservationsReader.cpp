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
 * MVIRIObservationsReader.cpp
 *
 *  Created on: 14 Oct 2015
 *      Author: alix
 */

#include "MFGObservationsReader.h"
#include "../../Utils/Math.h"

using std::string;
using std::vector;
using std::shared_ptr;
using std::make_shared;

using netCDF::NcVar;
using netCDF::NcFile;
using netCDF::exceptions::NcException;

using GEDAP::sqr;

/* ****************************************************************************************************	*
 * 											PUBLIC METHODS												*
 * ****************************************************************************************************	*/

MFGObservationsReader::MFGObservationsReader(shared_ptr<AbstractLogger> logger)
        : ObservationsReader(logger),
          sza(IMAGE_NB_LINES, IMAGE_NB_COLUMNS),
          saa(IMAGE_NB_LINES, IMAGE_NB_COLUMNS),
          vza(IMAGE_NB_LINES, IMAGE_NB_COLUMNS),
          vaa(IMAGE_NB_LINES, IMAGE_NB_COLUMNS),
          dq_bitmask(IMAGE_NB_LINES, IMAGE_NB_COLUMNS),
          acq_time(MAX_NB_BANDS),
          BRF(MAX_NB_BANDS),
          radiometric_noise(MAX_NB_BANDS){

    // Allocate memory for per-channel data
    for (auto &it : acq_time)
        it.resize(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
    for (auto &it : BRF)
        it.resize(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
    for (auto &it : radiometric_noise)
        it.resize(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);

}

void MFGObservationsReader::setup(const string imageFilename, shared_ptr<AbstractLogger> logger) {

    // Try to open input files
    isReady = true;

    try {
        logger->write_log(LogLevel::_DEBUG_, "MFGObservationsReader", "setup",
                          "Opening NetCDF image file " + imageFilename);
        mp_imageFile = make_shared<NcFile>(imageFilename, NcFile::read);

        vector<string> required_fields = {"solar_zenith_angle", "solar_azimuth_angle",
                                          "satellite_zenith_angle", "satellite_azimuth_angle",
                                          "sub_satellite_latitude_start", "sub_satellite_latitude_end",
                                          "sub_satellite_longitude_start", "sub_satellite_longitude_end",
                                          "u_independent_toa_bidirectional_reflectance","u_structured_toa_bidirectional_reflectance",
                                          "u_common_toa_bidirectional_reflectance",
                                          "toa_bidirectional_reflectance_vis", "time"};

        for (const string &field_name : required_fields) {
            if (mp_imageFile->getVar(field_name).isNull()) {
                logger->write_log(LogLevel::_ERROR_, "MFGObservationsReader", "setup",
                                  string("Field ").append(field_name).append(" is missing from image file")
                                  .append(imageFilename).append("; tile generation is impossible"));
                isReady = false;
            }
        }
    }
    catch (NcException &e) {
        /* NetCDF-CXX 4.3 required
        logger->write_log(LogLevel::_ERROR_, "MFGObservationsReader", "setup",
                          "Error while opening NetCDF file " + imageFilename +
                          " [NetCDF err code " + to_string(e.errorCode()) + "]");*/
        logger->write_log(LogLevel::_ERROR_, "MFGObservationsReader", "setup",
                          "Error while opening NetCDF file " + imageFilename);
        isReady = false;
    }
    catch (...) {
        logger->write_log(LogLevel::_ERROR_, "MFGObservationsReader", "setup",
                          "Error while opening " + imageFilename + "NetCDF file. Image skipped.");
        isReady = false;
    }

}


void MFGObservationsReader::read() {

    // Set channel to 0 (because there is only 1 channel in the case of MVIRI)
    // TODO: remove multi-channel facilities
    channel = 0;

    // Read acquisition time
    read_acquisition_time();

    // Read acquisition time
    read_data_quality_bitmask();

    // Read angles
    read_angles();

    // Read TOA BRF
    read_toa_brf();

    // Read error data
    read_error();

}


void MFGObservationsReader::extract_data(Observations &obs_data, TileDimensions &tile_dims) {
    // Extract angles
    extract_angles(obs_data, tile_dims);

    // Read rest of data (varies vs channel)
    for (channel = 0; channel < v_channels.size(); channel++) {
        extract_channel_observations(obs_data.v_band_obs[channel], tile_dims);
    }
}

/* ****************************************************************************************************	*
 * 											PRIVATE METHODS												*
 * ****************************************************************************************************	*/

void MFGObservationsReader::read_acquisition_time() {

    double add_offset = 0.0;

    // reset time array
    acq_time[channel] = Eigen::ArrayXXf_RowMaj::Zero(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);

    // Allocate buffer: the time variable is given at a different resolution (x_ir, y_ir)
    // TODO: add automatic buffer resizing and type detection
    Eigen::ArrayXXui_RowMaj packed_values(IMAGE_NB_LINES / 2, IMAGE_NB_COLUMNS / 2);

    try {
        // Read time data
        NcVar ncTime = mp_imageFile->getVar("time");
        ncTime.getVar(packed_values.data());
        ncTime.getAtt("add_offset").getValues(&add_offset);

        // Unpack
        for (int line = 0; line < IMAGE_NB_LINES; line++) {
            for (int column = 0; column < IMAGE_NB_COLUMNS; column++) {
                acq_time[channel](line, column) =
                        40587.0 + //days between November 17th 1858 and January 1st 1970
                        (packed_values(int(floor(line/ 2)), int(floor(column / 2))) + add_offset) /
                        (60 * 60 * 24);
            }
        }

        acq_time_info.min = acq_time[channel].minCoeff();
        acq_time_info.max = acq_time[channel].maxCoeff();
        acq_time_info.mean = acq_time[channel].mean();
        /*
        logger->write_log(LogLevel::_DEBUG_, "MFGObservationsReader", "read_acquisition_time",
                          "acq_time_info.min = " + to_string(acq_time_info.min));
        logger->write_log(LogLevel::_DEBUG_, "MFGObservationsReader", "read_acquisition_time",
                          "acq_time_info.max = " + to_string(acq_time_info.max));
        logger->write_log(LogLevel::_DEBUG_, "MFGObservationsReader", "read_acquisition_time",
                          "acq_time_info.mean = " + to_string(acq_time_info.mean));
          */



    }
    catch (NcException &c) {
        /* NetCDF-CXX 4.3 required
        logger->write_log(LogLevel::_ERROR_, "MFGObservationsReader", "read_acquisition_time",
                          "Error when reading NetCDF image file [NetCDF err code " +
                          std::to_string(c.errorCode()) + "]");*/
        logger->write_log(LogLevel::_ERROR_, "MFGObservationsReader", "read_acquisition_time",
                          "Error when reading NetCDF image file [NetCDF err code "); // FIXME: not very helpful: what file?
    }
}


void MFGObservationsReader::read_angles() {

    // Reset angle arrays
    sza = Eigen::ArrayXXf_RowMaj::Zero(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
    saa = Eigen::ArrayXXf_RowMaj::Zero(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
    vaa = Eigen::ArrayXXf_RowMaj::Zero(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
    vza = Eigen::ArrayXXf_RowMaj::Zero(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);

    // Read sun angles
    try {
        // Allocate buffers
        // the sun zenith and azimuth angles and the satellite zenith and azimuth angles are given at
        // a lower resolution x_tie, y_tie = (500, 500)

        // NOTE: Write a function to avoid rewriting the same code multiple times

        Eigen::ArrayXXus_RowMaj packed_values_ushort(IMAGE_NB_LINES/10, IMAGE_NB_COLUMNS/10);
        Eigen::ArrayXXs_RowMaj packed_values_short(IMAGE_NB_LINES/10, IMAGE_NB_COLUMNS/10);

        double scale_factor = 0.0;
        double add_offset = 0.0;
        double fill_value = 0.0;

        int x_1,x_2,y_1,y_2;


        // Read values for SZA
        NcVar ncSZA = mp_imageFile->getVar("solar_zenith_angle");
        ncSZA.getAtt("add_offset").getValues(&add_offset);
        ncSZA.getAtt("scale_factor").getValues(&scale_factor);
        ncSZA.getAtt("_FillValue").getValues(&fill_value);
        ncSZA.getVar(packed_values_short.data());


        // Unpack
        for (int line = 0; line < IMAGE_NB_LINES-1; line++){
            for (int column = 0; column < IMAGE_NB_COLUMNS-1; column++){
                if (line < 9 || column < 9 ){
                    sza(line, column) = -999.0;
                }
                else if ( line > 4980 || column > 4980){
                    sza(line, column) =  -999.0;
                }

                else {
                    // bilinear interpolation with
                    std::div_t dv_l{};
                    std::div_t dv_c{};

                    dv_l = std::div(line, 10);
                    dv_c = std::div(column, 10);

                    // x_1, x_2, y_1, y_2 are in the IMAGE_NB_LINES X IMAGE_NB_COLUMNS space (5000 X 5000)
                    // while x_tie and y_tie are in the IMAGE_NB_LINES/10 space (500 X 500)

                    int y_tie = dv_l.quot;
                    int x_tie = dv_c.quot;

                    x_1= (dv_c.quot*10)-1;
                    x_2= ((dv_c.quot+1)*10)-1;
                    y_1= (dv_l.quot*10)-1 ;
                    y_2= ((dv_l.quot+1)*10)-1;

                    float diff_x2_x1 = x_2-x_1;
                    float diff_y2_y1 = y_2-y_1;
                    float x2_column = x_2-column;
                    float column_x1 = column-x_1;

                    float k1 = ((x_2-column)/diff_x2_x1);
                    float k2 = ((column-x_1)/diff_x2_x1);

                    float Q11 = packed_values_short(y_tie,x_tie);
                    float Q21 = packed_values_short(y_tie,x_tie+1);
                    float Q12 = packed_values_short(y_tie+1,x_tie);
                    float Q22 = packed_values_short(y_tie+1,x_tie+1);

                    if ((Q11 != fill_value) && (Q21!=fill_value) && (Q12!=fill_value) && (Q22!=fill_value)){

                        Q11 = Q11*scale_factor + add_offset;
                        Q21 = Q21*scale_factor + add_offset;
                        Q12 = Q12*scale_factor + add_offset;
                        Q22 = Q22*scale_factor + add_offset;


                        float f_x_y1 = k1*(Q11) + k2*(Q21) ;
                        float f_x_y2 = k1*(Q12) + k2*(Q22);

                        sza(line, column) = (((y_2-line)/diff_y2_y1)*f_x_y1)+(((line-y_1)/diff_y2_y1)*f_x_y2);

                    }
                    else {
                        sza(line, column) = fill_value;
                    }
                }
            }
        }



        // Read values for SAA
        NcVar ncSAA = mp_imageFile->getVar("solar_azimuth_angle"); // Unsigned
        ncSAA.getAtt("add_offset").getValues(&add_offset);
        ncSAA.getAtt("scale_factor").getValues(&scale_factor);
        ncSAA.getAtt("_FillValue").getValues(&fill_value);
        ncSAA.getVar(packed_values_ushort.data());

        // Unpack
        for (int line = 0; line < IMAGE_NB_LINES-1; line++){
            for (int column = 0; column < IMAGE_NB_COLUMNS-1; column++){
                if (line < 9 || column < 9 ){
                    saa(line, column) = -999.0;
                }
                else if ( line > 4980 || column > 4980){
                    saa(line, column) =  -999.0;
                }

                else {
                    // bilinear interpolation with
                    std::div_t dv_l{};
                    std::div_t dv_c{};

                    dv_l = std::div(line, 10);
                    dv_c = std::div(column, 10);

                    // x_1, x_2, y_1, y_2 are in the IMAGE_NB_LINES X IMAGE_NB_COLUMNS space (5000 X 5000)
                    // while x_tie and y_tie are in the IMAGE_NB_LINES/10 space (500 X 500)

                    int y_tie = dv_l.quot;
                    int x_tie = dv_c.quot;

                    x_1= (dv_c.quot*10)-1;
                    x_2= ((dv_c.quot+1)*10)-1;
                    y_1= (dv_l.quot*10)-1 ;
                    y_2= ((dv_l.quot+1)*10)-1;

                    float diff_x2_x1 = x_2-x_1;
                    float diff_y2_y1 = y_2-y_1;
                    float x2_column = x_2-column;
                    float column_x1 = column-x_1;

                    float k1 = ((x_2-column)/diff_x2_x1);
                    float k2 = ((column-x_1)/diff_x2_x1);

                    float Q11 = packed_values_ushort(y_tie,x_tie);
                    float Q21 = packed_values_ushort(y_tie,x_tie+1);
                    float Q12 = packed_values_ushort(y_tie+1,x_tie);
                    float Q22 = packed_values_ushort(y_tie+1,x_tie+1);

                    if ((Q11 != fill_value) && (Q21!=fill_value) && (Q12!=fill_value) && (Q22!=fill_value)){

                        Q11 = Q11*scale_factor + add_offset;
                        Q21 = Q21*scale_factor + add_offset;
                        Q12 = Q12*scale_factor + add_offset;
                        Q22 = Q22*scale_factor + add_offset;


                        float f_x_y1 = k1*(Q11) + k2*(Q21) ;
                        float f_x_y2 = k1*(Q12) + k2*(Q22);

                        saa(line, column) = (((y_2-line)/diff_y2_y1)*f_x_y1)+(((line-y_1)/diff_y2_y1)*f_x_y2);

                        if ((Q11>300.0) && (Q12>300.0) && (Q21<50.0) && (Q22<50.0)) {
                            Q11 = 0;
                            Q12 = 0;

                            float f_x_y1 = k1*(Q11) + k2*(Q21) ;
                            float f_x_y2 = k1*(Q12) + k2*(Q22);

                            saa(line, column) = (((y_2-line)/diff_y2_y1)*f_x_y1)+(((line-y_1)/diff_y2_y1)*f_x_y2);
                        }
                        if ((Q11>300.0) && (Q12>300.0) && (Q21>300.0) && (Q22<50.0)) {
                            Q22 = 360.0;

                            float f_x_y1 = k1*(Q11) + k2*(Q21) ;
                            float f_x_y2 = k1*(Q12) + k2*(Q22);

                            saa(line, column) = (((y_2-line)/diff_y2_y1)*f_x_y1)+(((line-y_1)/diff_y2_y1)*f_x_y2);
                        }
                        if ((Q11>300.0) && (Q12<50.0) && (Q21<50.0) && (Q22<50.0)) {
                            Q11 = 0;

                            float f_x_y1 = k1*(Q11) + k2*(Q21) ;
                            float f_x_y2 = k1*(Q12) + k2*(Q22);

                            saa(line, column) = (((y_2-line)/diff_y2_y1)*f_x_y1)+(((line-y_1)/diff_y2_y1)*f_x_y2);
                        }


                    }
                    else {
                        saa(line, column) = fill_value;
                    }
                }
            }
        }



        Eigen::ArrayXXus_RowMaj packed_values_vza(IMAGE_NB_LINES/10, IMAGE_NB_COLUMNS/10);
        Eigen::ArrayXXus_RowMaj packed_values_vaa(IMAGE_NB_LINES/10, IMAGE_NB_COLUMNS/10);

        // Read values for VZA
        NcVar ncVZA = mp_imageFile->getVar("satellite_zenith_angle");
        ncVZA.getAtt("add_offset").getValues(&add_offset);
        ncVZA.getAtt("scale_factor").getValues(&scale_factor);
        ncVZA.getAtt("_FillValue").getValues(&fill_value);
        ncVZA.getVar(packed_values_vza.data());

        // Unpack
        for (int line = 0; line < IMAGE_NB_LINES-1; line++){
            for (int column = 0; column < IMAGE_NB_COLUMNS-1; column++){
                if (line < 9 || column < 9 ){
                    vza(line, column) = -999.0;
                }
                else if ( line > 4980 || column > 4980){
                    vza(line, column) =  -999.0;
                }

                else {
                    // bilinear interpolation with
                    std::div_t dv_l{};
                    std::div_t dv_c{};

                    dv_l = std::div(line, 10);
                    dv_c = std::div(column, 10);

                    // x_1, x_2, y_1, y_2 are in the IMAGE_NB_LINES X IMAGE_NB_COLUMNS space (5000 X 5000)
                    // while x_tie and y_tie are in the IMAGE_NB_LINES/10 space (500 X 500)

                    int y_tie = dv_l.quot;
                    int x_tie = dv_c.quot;

                    x_1= (dv_c.quot*10)-1;
                    x_2= ((dv_c.quot+1)*10)-1;
                    y_1= (dv_l.quot*10)-1 ;
                    y_2= ((dv_l.quot+1)*10)-1;

                    float diff_x2_x1 = x_2-x_1;
                    float diff_y2_y1 = y_2-y_1;
                    float x2_column = x_2-column;
                    float column_x1 = column-x_1;

                    float k1 = ((x_2-column)/diff_x2_x1);
                    float k2 = ((column-x_1)/diff_x2_x1);

                    float Q11 = packed_values_vza(y_tie,x_tie);
                    float Q21 = packed_values_vza(y_tie,x_tie+1);
                    float Q12 = packed_values_vza(y_tie+1,x_tie);
                    float Q22 = packed_values_vza(y_tie+1,x_tie+1);

                    if ((Q11 != fill_value) && (Q21!=fill_value) && (Q12!=fill_value) && (Q22!=fill_value)){

                        Q11 = Q11*scale_factor + add_offset;
                        Q21 = Q21*scale_factor + add_offset;
                        Q12 = Q12*scale_factor + add_offset;
                        Q22 = Q22*scale_factor + add_offset;


                        float f_x_y1 = k1*(Q11) + k2*(Q21) ;
                        float f_x_y2 = k1*(Q12) + k2*(Q22);

                        vza(line, column) = (((y_2-line)/diff_y2_y1)*f_x_y1)+(((line-y_1)/diff_y2_y1)*f_x_y2);
                    }
                    else {
                        vza(line, column) = fill_value;
                    }
                }
            }
        }


        // Read values for VAA
        NcVar ncVAA = mp_imageFile->getVar("satellite_azimuth_angle"); // Unsigned
        ncVAA.getAtt("add_offset").getValues(&add_offset);
        ncVAA.getAtt("scale_factor").getValues(&scale_factor);
        ncVAA.getAtt("_FillValue").getValues(&fill_value);
        ncVAA.getVar(packed_values_vaa.data());

        // Unpack
        for (int line = 0; line < IMAGE_NB_LINES-1; line++){
            for (int column = 0; column < IMAGE_NB_COLUMNS-1; column++){
                if (line < 9 || column < 9 ){
                    vaa(line, column) = -999.0;
                }
                else if ( line > 4980 || column > 4980){
                    vaa(line, column) =  -999.0;
                }

                else {
                    // bilinear interpolation with
                    std::div_t dv_l{};
                    std::div_t dv_c{};

                    dv_l = std::div(line, 10);
                    dv_c = std::div(column, 10);

                    // x_1, x_2, y_1, y_2 are in the IMAGE_NB_LINES X IMAGE_NB_COLUMNS space (5000 X 5000)
                    // while x_tie and y_tie are in the IMAGE_NB_LINES/10 space (500 X 500)

                    int y_tie = dv_l.quot;
                    int x_tie = dv_c.quot;

                    x_1= (dv_c.quot*10)-1;
                    x_2= ((dv_c.quot+1)*10)-1;
                    y_1= (dv_l.quot*10)-1 ;
                    y_2= ((dv_l.quot+1)*10)-1;

                    float diff_x2_x1 = x_2-x_1;
                    float diff_y2_y1 = y_2-y_1;
                    float x2_column = x_2-column;
                    float column_x1 = column-x_1;

                    float k1 = ((x_2-column)/diff_x2_x1);
                    float k2 = ((column-x_1)/diff_x2_x1);

                    float Q11 = packed_values_vaa(y_tie,x_tie);
                    float Q21 = packed_values_vaa(y_tie,x_tie+1);
                    float Q12 = packed_values_vaa(y_tie+1,x_tie);
                    float Q22 = packed_values_vaa(y_tie+1,x_tie+1);

                    if ((Q11 != fill_value) && (Q21!=fill_value) && (Q12!=fill_value) && (Q22!=fill_value)){

                        Q11 = Q11*scale_factor + add_offset;
                        Q21 = Q21*scale_factor + add_offset;
                        Q12 = Q12*scale_factor + add_offset;
                        Q22 = Q22*scale_factor + add_offset;


                        float f_x_y1 = k1*(Q11) + k2*(Q21) ;
                        float f_x_y2 = k1*(Q12) + k2*(Q22);

                        vaa(line, column) = (((y_2-line)/diff_y2_y1)*f_x_y1)+(((line-y_1)/diff_y2_y1)*f_x_y2);

                        if ((Q11>300.0) && (Q12>300.0) && (Q21<50.0) && (Q22<50.0)) {
                            Q11 = 0;
                            Q12 = 0;

                            float f_x_y1 = k1*(Q11) + k2*(Q21) ;
                            float f_x_y2 = k1*(Q12) + k2*(Q22);

                            vaa(line, column) = (((y_2-line)/diff_y2_y1)*f_x_y1)+(((line-y_1)/diff_y2_y1)*f_x_y2);
                        }
                    }
                    else {
                        vaa(line, column) = fill_value;
                    }
                }
            }
        }
    }

    catch (NcException &c) {
        /* NetCDF-CXX 4.3 required
        logger->write_log(LogLevel::_ERROR_, "MFGObservationsReader", "read_angles",
                          "Error when reading NetCDF image file [NetCDF err code " +
                          std::to_string(c.errorCode()) + "]");*/
        logger->write_log(LogLevel::_ERROR_, "MFGObservationsReader", "read_angles",
                          "Error when reading NetCDF image file");
    }
}


void MFGObservationsReader::read_toa_brf() {
    double scale_factor = 0.0;
    double add_offset = 0.0;

    // Reset BRF array
    BRF[channel] = Eigen::ArrayXXf_RowMaj::Zero(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);

    // Allocate buffer
    Eigen::ArrayXXus_RowMaj packed_values(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);

    try {
        // Read TOA BRF data
        NcVar ncTOA = mp_imageFile->getVar("toa_bidirectional_reflectance_vis");
        ncTOA.getVar(packed_values.data());
        ncTOA.getAtt("add_offset").getValues(&add_offset);
        ncTOA.getAtt("scale_factor").getValues(&scale_factor);

        // Unpack
        for (int line = 0; line < IMAGE_NB_LINES; line++) {
            for (int column = 0; column < IMAGE_NB_COLUMNS; column++) {
                BRF[channel](line, column) = (packed_values(line, column) * scale_factor) + add_offset;
            }
        }
    }
    catch (NcException &c) {
        /* NetCDF 4.3 required
        logger->write_log(LogLevel::_ERROR_,
                          "MFGObservationsReader", "read_toa_brf",
                          "Error when reading NetCDF image file [NetCDF err code " +
                          std::to_string(c.errorCode()) + "]");*/
        logger->write_log(LogLevel::_ERROR_,
                          "MFGObservationsReader", "read_toa_brf",
                          "Error when reading MVIRI image NetCDF file"); //
    }

}


void MFGObservationsReader::read_error() {

    // Reset radiometric noise array
    radiometric_noise[channel] = Eigen::ArrayXXf_RowMaj::Zero(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);

    // Contribution due to uncertainty on BRF
    // --- Allocate buffers and temporaries
    Eigen::ArrayXXus_RowMaj packed_values_random(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
    Eigen::ArrayXXus_RowMaj packed_values(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
    double add_offset_random = 0., scale_factor_random = 0., add_offset = 0., scale_factor = 0., sigma_i_squared = 0.;

    // --- Load image data
    try {
        // Get random uncertainty
        mp_imageFile->getVar("u_independent_toa_bidirectional_reflectance").getVar(packed_values_random.data());

        // Get values of scale factor and offset
        mp_imageFile->getVar("u_independent_toa_bidirectional_reflectance").getAtt("add_offset").getValues(&add_offset_random);
        mp_imageFile->getVar("u_independent_toa_bidirectional_reflectance").getAtt("scale_factor").getValues(&scale_factor_random);

        // Get non random uncertainty
        mp_imageFile->getVar("u_structured_toa_bidirectional_reflectance").getVar(packed_values.data());

        // Get values of scale factor and offset
        mp_imageFile->getVar("u_structured_toa_bidirectional_reflectance").getAtt("add_offset").getValues(&add_offset);
        mp_imageFile->getVar("u_structured_toa_bidirectional_reflectance").getAtt("scale_factor").getValues(&scale_factor);

//        // Get common uncertainty
//        mp_imageFile->getVar("u_common_toa_bidirectional_reflectance").getVar(packed_values.data());

//        // Get values of scale factor and offset
//        mp_imageFile->getVar("u_common_toa_bidirectional_reflectance").getAtt("add_offset").getValues(&add_offset);
//        mp_imageFile->getVar("u_common_toa_bidirectional_reflectance").getAtt("scale_factor").getValues(&scale_factor);

    }
    catch (NcException &c) {
        /* NetCDF-CXX 4.3 required
        logger->write_log(LogLevel::_ERROR_, "MFGObservationsReader", "read_error",
                          "Error when reading NetCDF image file [NetCDF err code " +
                          std::to_string(c.errorCode()) + "]");*/
        logger->write_log(LogLevel::_ERROR_, "MFGObservationsReader", "read_error",
                          "Error when reading NetCDF image file");
    }

    // --- Unpack data and apply multiplier to get absolute uncertainty
    for (int line = 0; line < IMAGE_NB_LINES; line++) {
        for (int column = 0; column < IMAGE_NB_COLUMNS; column++) {

            sigma_i_squared =
                    sqr((packed_values_random(line, column) * scale_factor_random + add_offset_random) * BRF[channel](line, column)) +
                    sqr((packed_values(line, column) * scale_factor + add_offset) * BRF[channel](line, column));
            radiometric_noise[channel](line, column) += sigma_i_squared;
        }
    }

    // Add contribution due to uncertainty on geolocation / coregistration
    // The derivatives are computed with a 2nd order centered finite difference approximation
    // For now, edges (where a backward / forward approximation is required) are omitted, and there is therefore a
    // 1-pixel margin around the image where the total uncertainty is wrong; a proper partial differentiating function
    // should be added to the math library
    // --- Allocate buffers and temporaries
    double sigma_x = 0., sigma_y = 0., sigma_r_squared = 0.;

    // --- Get uncertainty on position (per pixel)
    try {
        mp_imageFile->getVar("toa_bidirectional_reflectance_vis").getAtt("rms_landmarks_x_vis").getValues(&sigma_x);
        mp_imageFile->getVar("toa_bidirectional_reflectance_vis").getAtt("rms_landmarks_y_vis").getValues(&sigma_y);
    }
    catch (NcException &c) {
        /* NetCDF-CXX 4.3 required
        logger->write_log(LogLevel::_ERROR_, "MFGObservationsReader", "read_error",
                          "Error when reading NetCDF image file [NetCDF err code " +
                          std::to_string(c.errorCode()) + "]");*/
        logger->write_log(LogLevel::_ERROR_, "MFGObservationsReader", "read_error",
                          "Error when reading NetCDF image file");
    }

    // Computation of derivative of BRF: since sigma_x and sigma_y are given per pixel,
    // the central derivative requires h=2 (spacing between the extremes, in units of pixels)
    // lines correspond to y values while columns correspond to x values

    for (int line = 1; line < IMAGE_NB_LINES-1; ++line) {
        for (int column = 1; column < IMAGE_NB_COLUMNS-1; ++column) {

            sigma_r_squared =
                    sqr((BRF[channel](line + 1, column) - BRF[channel](line - 1, column)) / 2. * sigma_y) +
                    sqr((BRF[channel](line, column + 1) - BRF[channel](line, column - 1)) / 2. * sigma_x);
            radiometric_noise[channel](line, column) += sigma_r_squared;
        }
    }
}

void MFGObservationsReader::extract_angles(Observations &obs_data, TileDimensions &tile_dims) {
    int i_image_line = 0, i_image_column = 0;

    for (int i_line = 0; i_line < obs_data.nb_lines; i_line++) {
        for (int i_column = 0; i_column < obs_data.nb_columns; i_column++) {
            i_image_line = tile_dims.first_line + i_line;
            i_image_column = tile_dims.first_column + i_column;

            obs_data.SZA(i_line, i_column) = sza(i_image_line, i_image_column);
            obs_data.SAA(i_line, i_column) = saa(i_image_line, i_image_column);
            obs_data.VZA(i_line, i_column) = vza(i_image_line, i_image_column);
            obs_data.VAA(i_line, i_column) = vaa(i_image_line, i_image_column);

        }
    }
}

void MFGObservationsReader::extract_channel_observations(BandObservations &band_obs_data, TileDimensions &tile_dims) {
    int i_image_line = 0, i_image_column = 0;

    for (int i_line = 0; i_line < tile_dims.nb_lines; ++i_line) {
        for (int i_column = 0; i_column < tile_dims.nb_columns; ++i_column) {
            i_image_line = tile_dims.first_line + i_line;
            i_image_column = tile_dims.first_column + i_column;

            band_obs_data.acquisitionTime(i_line, i_column) =
                    acq_time[channel](i_image_line, i_image_column);
            band_obs_data.TOABRF(i_line, i_column) =
                    BRF[channel](i_image_line, i_image_column);
            // this value is passed to CISAR as sigma_TOA_BRF. Therefore it has to be a sigma and not a variance
            // the square root is mandatory here
            band_obs_data.radiometricNoise(i_line, i_column) =
                    sqrt(radiometric_noise[channel](i_image_line, i_image_column));

        }
    }
}

void MFGObservationsReader::read_data_quality_bitmask(){
    dq_bitmask = Eigen::ArrayXXi_RowMaj::Zero(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);

    Eigen::MatrixXi_RowMaj packed_values(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);

    // read and unpack data_quality_bitmask
    mp_imageFile->getVar("data_quality_bitmask").getVar(packed_values.data());
    for (int line = 0; line < IMAGE_NB_LINES; line++) {
        for (int column = 0; column < IMAGE_NB_COLUMNS; column++) {
            dq_bitmask(line, column) = packed_values(line,column);
        }
    }
}


void MFGObservationsReader::extract_data_quality_bitmask(Mask &mask_data, TileDimensions &tile_dims) {
    int i_image_line, i_image_column;
    for (int i_line = 0; i_line < mask_data.nb_lines; i_line++) {
        for (int i_column = 0; i_column < mask_data.nb_columns; i_column++) {
            i_image_line = tile_dims.first_line + i_line;
            i_image_column = tile_dims.first_column + i_column;

            bool isBad = false;
            isBad = (dq_bitmask(i_image_line, i_image_column) != 0);

            std::bitset<Mask::NB_MASK_VALUES> maskValue;

            if (isBad) maskValue.set(Mask::QUALITY);

            mask_data.data_quality_bitmask(i_line, i_column) = Mask::bitset_to_uchar(maskValue);
//            std::cout << "input_data.data_quality_bitmask(i_line, i_column)  in MFGObservationsReader =   " << std::bitset<6>(mask_data.data_quality_bitmask(i_line, i_column)) << std::endl;

        }
    }
}
