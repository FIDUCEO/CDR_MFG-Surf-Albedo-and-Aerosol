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
 * ModelParamsReader.cpp
 *
 *  Created on: 13 Oct 2015
 *      Author: alix
 */

#include "ModelParamsReader.h"

using namespace std;

using netCDF::NcFile;
using netCDF::exceptions::NcException;

void ModelParamsReader::extract_data(ModelParams &model_params_data, const string elevation_filename,
                                     const string geopotential_filename, const string aerosol_height_filename_1,
                                     const string aerosol_height_filename_2, const int dayOfMonth, const bool heightIsConstant,
                                     TileDimensions &tile_dims, const int step, const bool ecmwf_constant)
{
// these are the dimensions of the input tile we have to create
    int first_line   = tile_dims.first_line;
    int first_column = tile_dims.first_column;
    int nb_lines     = tile_dims.nb_lines;
    int nb_columns   = tile_dims.nb_columns;

    Eigen::MatrixXf_RowMaj U10_interp = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
    Eigen::MatrixXf_RowMaj V10_interp = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);

    Eigen::MatrixXf_RowMaj surfPressTemp = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
    Eigen::MatrixXf_RowMaj TCWVTemp = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
    Eigen::MatrixXf_RowMaj aerosolHeightTemp_1 = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
    Eigen::MatrixXf_RowMaj aerosolHeightTemp_2 = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
    Eigen::MatrixXf_RowMaj aerosol_elevation = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);


// model_params_data are the actual fields in the input tile that need to be filled
    if(ecmwf_constant){
        model_params_data.TCO3.setConstant(310.0);
        TCWVTemp.setConstant(15.0);
        surfPressTemp.setConstant(1013.0);
        model_params_data.windSpeed.setConstant(3.0);
        model_params_data.windDir.setConstant(0.0);
    }
    else {

// at this point we have a full ERA-INTERIM image (stepData) which has 241 lines (latitude) and 480 columns (longitude)
// and has been correctly time interpolated. Now we need to perform the pixel spatial interpolation => we need to fill the input tile data
// with the correct interpolated info.

        // Update cache if necessary
        m_geopotential_filename = geopotential_filename;
        m_elevation_filename = elevation_filename;
        m_aerosol_height_filename_1 = aerosol_height_filename_1;
        m_aerosol_height_filename_2 = aerosol_height_filename_2;

        this->update_cache(ecmwf_constant);

//        pick one generic pixel to evaluate the half width between pixels
        int i=20;
        float half_width= (m_era_interim_longitude_cache(i+1)-m_era_interim_longitude_cache(i))/2;

        vector<int> latitude_indices;
        vector<int> longitude_indices;
        int i_pixel = 0;
        for(int i_line = 0; i_line < nb_lines; i_line++)
        {
            for(int i_column = 0; i_column < nb_columns; i_column++)
            {
                if (m_latitude_cache(first_line+i_line, first_column+i_column) != -999 && m_longitude_cache(first_line+i_line, first_column+i_column) != -999){
                    int latIndex = 0;
                    int lonIndex = 0;

                    while(m_era_interim_latitude_cache(latIndex) > m_latitude_cache(first_line+i_line, first_column+i_column))
                    {
                        latIndex ++;
                    }
                    while(m_era_interim_longitude_cache(lonIndex) < m_longitude_cache(first_line+i_line, first_column+i_column))
                    {
                        lonIndex ++;
                    }

                    latIndex = latIndex-1;
                    lonIndex = lonIndex-1;

                    latitude_indices.push_back(latIndex);
                    longitude_indices.push_back(lonIndex);

                    if (m_latitude_cache(first_line+i_line, first_column+i_column) > (m_era_interim_latitude_cache(latIndex)-half_width))
                    {
                        if (latIndex==0){latitude_indices.push_back(latIndex);}
                        else {latitude_indices.push_back(latIndex-1);}
                    }
                    else{
                        if (latIndex==ERA_INTERIM_IMAGE_NB_LINES) {latitude_indices.push_back(latIndex);}
                        else {latitude_indices.push_back(latIndex+1);}
                    }

                    if (m_longitude_cache(first_line+i_line, first_column+i_column) < (m_era_interim_longitude_cache(latIndex)+half_width))
                    {
                        if (lonIndex==0) {longitude_indices.push_back(lonIndex);}
                        else {longitude_indices.push_back(lonIndex-1);}
                    }
                    else{
                        if (lonIndex==ERA_INTERIM_IMAGE_NB_COLUMNS) {longitude_indices.push_back(lonIndex); }
                        else {longitude_indices.push_back(lonIndex+1);}

                    }

                    float x_1 = m_era_interim_longitude_cache(longitude_indices.at(0));
                    float x_2 = m_era_interim_longitude_cache(longitude_indices.at(1));

                    float y_1 = m_era_interim_latitude_cache(latitude_indices.at(0));
                    float y_2 = m_era_interim_latitude_cache(latitude_indices.at(1));

                    float line = m_latitude_cache(first_line+i_line, first_column+i_column);
                    float column = m_longitude_cache(first_line+i_line, first_column+i_column);

                    float diff_x2_x1 = x_2 - x_1;
                    float diff_y2_y1 = y_2 - y_1;

                    float x2_column = x_2 - column;
                    float column_x1 = column - x_1;

                    float y2_line = y_2 - line;
                    float line_y1 = line - y_1;

                    float k1 = (x2_column/diff_x2_x1);
                    float k2 = (column_x1/diff_x2_x1);


                    model_params_data.TCO3(i_pixel) = (y2_line/diff_y2_y1) * (k1*TCO3.stepData(latitude_indices.at(0), longitude_indices.at(0)) + k2*TCO3.stepData(latitude_indices.at(0), longitude_indices.at(1)))   +
                            (line_y1/diff_y2_y1) * (k1*TCO3.stepData(latitude_indices.at(1), longitude_indices.at(0)) + k2*TCO3.stepData(latitude_indices.at(1), longitude_indices.at(1)));

                    TCWVTemp(i_pixel) =  (y2_line/diff_y2_y1) * (k1*TCWV.stepData(latitude_indices.at(0), longitude_indices.at(0)) + k2*TCWV.stepData(latitude_indices.at(0), longitude_indices.at(1)))  +
                            (line_y1/diff_y2_y1) * (k1*TCWV.stepData(latitude_indices.at(1), longitude_indices.at(0)) + k2*TCWV.stepData(latitude_indices.at(1), longitude_indices.at(1)));

                    surfPressTemp(i_pixel) = (y2_line/diff_y2_y1) * (k1*SrfPress.stepData(latitude_indices.at(0), longitude_indices.at(0)) + k2*SrfPress.stepData(latitude_indices.at(0), longitude_indices.at(1)))  +
                            (line_y1/diff_y2_y1) * (k1*SrfPress.stepData(latitude_indices.at(1), longitude_indices.at(0)) + k2*SrfPress.stepData(latitude_indices.at(1), longitude_indices.at(1)));

                    U10_interp(i_pixel) = (y2_line/diff_y2_y1) * (k1*U10_era_image.stepData(latitude_indices.at(0), longitude_indices.at(0)) + k2*U10_era_image.stepData(latitude_indices.at(0), longitude_indices.at(1)))  +
                            (line_y1/diff_y2_y1) * (k1*U10_era_image.stepData(latitude_indices.at(1), longitude_indices.at(0)) + k2*U10_era_image.stepData(latitude_indices.at(1), longitude_indices.at(1)));

                    V10_interp(i_pixel) = (y2_line/diff_y2_y1) * (k1*V10_era_image.stepData(latitude_indices.at(0), longitude_indices.at(0)) + k2*V10_era_image.stepData(latitude_indices.at(0), longitude_indices.at(1)))  +
                            (line_y1/diff_y2_y1) * (k1*V10_era_image.stepData(latitude_indices.at(1), longitude_indices.at(0)) + k2*V10_era_image.stepData(latitude_indices.at(1), longitude_indices.at(1)));

                    /* calculate wind speed and direction
                    * wind speed   	--> sqrt( u10^2 + v10^2 )
                    * wind direction	--> 180/pi * atan2(u,v) + 180
                    * V component		--> -cos(dir * pi/180) * speed
                    * U component		--> -sin(dir * pi/180) * speed
                    */

                    model_params_data.windSpeed(i_pixel) = sqrt( U10_interp(i_pixel)  * U10_interp(i_pixel)  + V10_interp(i_pixel) * V10_interp(i_pixel) );
                    model_params_data.windDir(i_pixel) = COEFF_RAD2DEG * atan2(U10_interp(i_pixel), V10_interp(i_pixel)) + 180.0;

                    longitude_indices.clear();
                    latitude_indices.clear();

                }
            i_pixel++;
            }
        }

        int x;
        int altIndex;
        float elevationGeoPot;
        float elevationStatic;
        float TCWVSP;
        float actualheight;

        // Levels in km for the US76 standard profile
        float altStandardProfile[68] = {
            120, 115, 110, 105, 100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 47.5, 45, 42.5, 40, 37.5, 35, 32.5, 30, 27.5,
             25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1.9, 1.8, 1.7, 1.6,
             1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0};

//        // Fraction of TCWV in each layer for a  US76 standard profile
        float StandardProfileTCWVFrac[68] = {
     0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000001,
     0.0000002, 0.0000004, 0.0000009, 0.0000018, 0.0000025, 0.0000034, 0.0000047, 0.0000066, 0.0000093, 0.0000131,
     0.0000185, 0.0000264, 0.0000378, 0.0000541, 0.0000623, 0.0000717, 0.0000823, 0.0000945, 0.0001085, 0.0001247,
     0.0001434, 0.0001653, 0.0001914, 0.0002267, 0.0002769, 0.0003687, 0.0005599, 0.0009735, 0.0018842, 0.0041090,
     0.0098794, 0.0213525, 0.0418681, 0.0773371, 0.1378394, 0.2386593, 0.4020520, 0.4226675, 0.4441861, 0.4666079,
     0.4899327, 0.5141607, 0.5392918, 0.5653260, 0.5922633, 0.6201037, 0.6488473, 0.6786344, 0.7096055, 0.7417607,
     0.7750999, 0.8096232, 0.8453305, 0.8822218, 0.9202972, 0.9595566, 1.0000000, 1.0};

        for(int a = 0; a < nb_lines; a++)
        {
            for(int b = 0; b < nb_columns; b++)
            {
                //determine the actual pixel elevation
                elevationStatic = m_elevation_cache(first_line + a, first_column + b);
                //and the elevation at which the TCWV and surface pressure are determined
                elevationGeoPot = m_geopotential_cache(first_line + a, first_column + b) /9.80665;

                if(elevationGeoPot < 0)
                {
                    elevationGeoPot = 0;
                }

                //determine the aerosol layer height
                aerosolHeightTemp_1(a,b) = m_aerosol_height_cache_1(first_line + a, first_column + b)  - elevationStatic/1000.;
                aerosolHeightTemp_2(a,b) = m_aerosol_height_cache_2(first_line + a, first_column + b)  - elevationStatic/1000.;

                //Rescale the surace pressure at the pixel elevation
                actualheight = elevationStatic - elevationGeoPot; //m
                surfPressTemp(a,b) = surfPressTemp(a,b)*exp(-(9.80665*0.0289644/(8.31447*288.15))*actualheight);

                //Rescale the TCWV according to the pixel elevation
                // As pixels are usually not very elevated, we use a very simple search technic
                // starting from the surface to find the value of altIndex
                altIndex = 67;

                while((altStandardProfile[altIndex] <= actualheight/1000.) && (altIndex > 0)) {
                    altIndex--;
                }

                // Rescale TCWV
                TCWVTemp(a,b) = TCWVTemp(a,b) * StandardProfileTCWVFrac[altIndex];
            }
        }

        model_params_data.surfPress = surfPressTemp;
        model_params_data.TCWV = TCWVTemp;

        // time interpolation of aerosol layer height. Depending on the day of the month, if day > 10 File_1 is the current month and File_2 is the next month
        // if day < 10 File_1 is the precedent month, Fi
        if(!heightIsConstant)
        {
            int x;
            if(dayOfMonth < 10)
            {
                x = 20 + dayOfMonth;
            }
            else
            {
                x = dayOfMonth - 10;
            }
            for(int a = 0; a < nb_lines; a++)
            {
                for(int b = 0; b < nb_columns; b++)
                {
                    aerosol_elevation(a,b) = (aerosolHeightTemp_2(a, b) - aerosolHeightTemp_1(a,b)) * x / 30 + aerosolHeightTemp_1(a,b);
                }
            }
            model_params_data.aerosolHeight = aerosol_elevation;
        }
        else
        {
            model_params_data.aerosolHeight.setConstant(2.0);
        }

    }
}

void ModelParamsReader::read_first_model_params_data(const string inputFile, const int ecmwf_image_slot)
{

    mp_inputFile = make_shared<NcFile>(inputFile, NcFile::read);

    m_era_interim_latitude_cache.resize(ERA_INTERIM_IMAGE_NB_LINES);
    m_era_interim_longitude_cache.resize(ERA_INTERIM_IMAGE_NB_COLUMNS);

    mp_inputFile->getVar("latitude").getVar(m_era_interim_latitude_cache.data());
    mp_inputFile->getVar("longitude").getVar(m_era_interim_longitude_cache.data());

    // read data for the whole image and save them into nextData array
    this->load_era_interim_cache(inputFile, ecmwf_image_slot);
}

void ModelParamsReader::compute_gradient_for_interpolation(const string inputFile, const int ecmwf_image_slot)
{
    // reset gradients
    TCO3.reset_gradient();
    TCWV.reset_gradient();
    SrfPress.reset_gradient();
    U10_era_image.reset_gradient();
    V10_era_image.reset_gradient();

    TCO3.set_next_data_to_current();
    TCWV.set_next_data_to_current();
    SrfPress.set_next_data_to_current();
    U10_era_image.set_next_data_to_current();
    V10_era_image.set_next_data_to_current();

    //Load the subsequent frame of the ERA-INTERIM image
    this->load_era_interim_cache(inputFile, ecmwf_image_slot);

    // compute and set gradients
    TCO3.compute_and_set_gradient();
    TCWV.compute_and_set_gradient();
    SrfPress.compute_and_set_gradient();
    U10_era_image.compute_and_set_gradient();
    V10_era_image.compute_and_set_gradient();
}

void ModelParamsReader::create_full_era_interim_image_at_step(const int slot)
{
    TCO3.evaluate_data_at_step(slot%12);
    TCWV.evaluate_data_at_step(slot%12);
    SrfPress.evaluate_data_at_step(slot%12);
    U10_era_image.evaluate_data_at_step(slot%12);
    V10_era_image.evaluate_data_at_step(slot%12);
}


/* ****************************************************************************************************	*
 * 											PRIVATE METHODS												*
 * ****************************************************************************************************	*/
void ModelParamsReader::load_era_interim_cache(const string inputFile, const int slot)
{
    double scale_factor = 0.0;
    double add_offset = 0.0;
    short missing_value = -32767;
    double scale_factor_bis = 0.0;
    double add_offset_bis = 0.0;

    try
    {
        mp_inputFile = make_shared<NcFile>(inputFile, NcFile::read);

        // Allocate buffers
        std::vector<size_t> start = {static_cast<size_t>(slot),
                                     0,
                                     0};
        std::vector<size_t> count = {1,
                                     static_cast<size_t>(ERA_INTERIM_IMAGE_NB_LINES),
                                     static_cast<size_t>(ERA_INTERIM_IMAGE_NB_COLUMNS)};


        packed_values.resize(ERA_INTERIM_IMAGE_NB_LINES, ERA_INTERIM_IMAGE_NB_COLUMNS);
        packed_values_bis.resize(ERA_INTERIM_IMAGE_NB_LINES, ERA_INTERIM_IMAGE_NB_COLUMNS);


        /* ****************************	*
         * 				TCO3			*
         * ****************************	*/
        mp_inputFile->getVar("tco3").getVar(start, count, packed_values.data());
        mp_inputFile->getVar("tco3").getAtt("scale_factor").getValues(&scale_factor);
        mp_inputFile->getVar("tco3").getAtt("add_offset").getValues(&add_offset);

        for (int line = 0; line < ERA_INTERIM_IMAGE_NB_LINES; line++){
            for (int column = 0; column < ERA_INTERIM_IMAGE_NB_COLUMNS; column++){
                TCO3.nextData(line, column) = (packed_values(line, column) != missing_value) ?
                                              ((packed_values(line, column) * scale_factor) + add_offset) *
                                              COEFF_KGM22DU : 0.0;
            }
        }

        /* ****************************	*
         * 				TCWV			*
         * ****************************	*/

        mp_inputFile->getVar("tcwv").getVar(start, count, packed_values.data());
        mp_inputFile->getVar("tcwv").getAtt("scale_factor").getValues(&scale_factor);
        mp_inputFile->getVar("tcwv").getAtt("add_offset").getValues(&add_offset);

        for(int line = 0; line < ERA_INTERIM_IMAGE_NB_LINES; line++){
            for(int column = 0; column < ERA_INTERIM_IMAGE_NB_COLUMNS; column++){
                TCWV.nextData(line, column) = (packed_values(line, column) != missing_value) ?
                                              (packed_values(line, column) * scale_factor) + add_offset : 0.0;
            }
        }

        /* ****************************	*
         * 		SURFACE PRESSURE		*
         * ****************************	*/

        mp_inputFile->getVar("sp").getVar(start, count, packed_values.data());
        mp_inputFile->getVar("sp").getAtt("scale_factor").getValues(&scale_factor);
        mp_inputFile->getVar("sp").getAtt("add_offset").getValues(&add_offset);

        for(int line = 0; line < ERA_INTERIM_IMAGE_NB_LINES; line++){
            for(int column = 0; column < ERA_INTERIM_IMAGE_NB_COLUMNS; column++){
                SrfPress.nextData(line, column) = (packed_values(line, column) != missing_value) ?
                                                  ((packed_values(line, column) * scale_factor) + add_offset) / 100.0
                                                                                                   : 0.0;
            }
        }

        /* ****************************	*
         * 		WIND DIR & SPEED		*
         * ****************************	*/
        mp_inputFile->getVar("u10").getVar(start, count, packed_values.data());
        mp_inputFile->getVar("u10").getAtt("scale_factor").getValues(&scale_factor);
        mp_inputFile->getVar("u10").getAtt("add_offset").getValues(&add_offset);

        mp_inputFile->getVar("v10").getVar(start, count, packed_values_bis.data());
        mp_inputFile->getVar("v10").getAtt("scale_factor").getValues(&scale_factor_bis);
        mp_inputFile->getVar("v10").getAtt("add_offset").getValues(&add_offset_bis);

        for(int line = 0; line < ERA_INTERIM_IMAGE_NB_LINES; line++){
            for(int column = 0; column < ERA_INTERIM_IMAGE_NB_COLUMNS; column++){
                U10_era_image.nextData(line, column) =  (packed_values(line, column) != missing_value) ? (packed_values(line, column) * scale_factor) + add_offset : 0.0;
                V10_era_image.nextData(line, column) =  (packed_values_bis(line, column) != missing_value) ? (packed_values_bis(line, column) * scale_factor_bis) + add_offset_bis : 0.0;

//                WindSpeed.nextData(line, column) = sqrt( U10 * U10 + V10 * V10 );
//                WindDir.nextData(line, column)   = COEFF_RAD2DEG * atan2(U10, V10) + 180.0;
            }
        }
    }
    catch (NcException& c){
        logger->write_log(LogLevel::_ERROR_, "ModelParamsReader", "constructor", "can't open or read ERA-INTERIM file with path " + inputFile);
    }
}


void ModelParamsReader::update_cache(const bool ecmwf_constant)
{
    logger->write_log(LogLevel::_DEBUG_, "ModelParamsReader", "update_cache",
                      "Updating static data cache");

    if (!m_cache_ok)
    {
        // Open netCDF data files

        // -- Geopotential
        shared_ptr<NcFile> p_geopotential_file;
        try {
            logger->write_log(LogLevel::_DEBUG_, "ModelParamsReader", "update_cache",
                              "Opening NetCDF geopotential file " + m_geopotential_filename);
            p_geopotential_file = make_shared<NcFile>(m_geopotential_filename, NcFile::read);
        }
        catch (NcException &c) {
            logger->write_log(LogLevel::_ERROR_, "ModelParamsReader", "update_cache",
                              "Error while attempting to open geopotential file " + m_geopotential_filename);
            throw c;
        }

        // -- Elevation
        shared_ptr<NcFile> p_elevation_file;
        try {
            logger->write_log(LogLevel::_DEBUG_, "ModelParamsReader", "update_cache",
                              "Opening NetCDF elevation file" + m_elevation_filename);
            p_elevation_file = make_shared<NcFile>(m_elevation_filename, NcFile::read);
        }
        catch (NcException &c) {
            logger->write_log(LogLevel::_ERROR_, "ModelParamsReader", "update_cache",
                              "Error while attempting to open elevation file " + m_elevation_filename);
            throw c;
        }

        // -- Aerosol height current month
        shared_ptr<NcFile> p_aerosol_height_file_1;
        try {
            logger->write_log(LogLevel::_DEBUG_, "ModelParamsReader", "update_cache",
                              "Opening first NetCDF aerosol height file " + m_aerosol_height_filename_1);
            p_aerosol_height_file_1 = make_shared<NcFile>(m_aerosol_height_filename_1, NcFile::read);
        }
        catch (NcException &c) {
            logger->write_log(LogLevel::_ERROR_, "ModelParamsReader", "update_cache",
                              "Error while attempting to open aerosol height file " + m_aerosol_height_filename_1);
            throw c;
        }

        // -- Aerosol height next month
        shared_ptr<NcFile> p_aerosol_height_file_2;
        try {
            logger->write_log(LogLevel::_DEBUG_, "ModelParamsReader", "update_cache",
                              "Opening second NetCDF aerosol height file " + m_aerosol_height_filename_2);
            p_aerosol_height_file_2 = make_shared<NcFile>(m_aerosol_height_filename_2, NcFile::read);
        }
        catch (NcException &c) {
            logger->write_log(LogLevel::_ERROR_, "ModelParamsReader", "update_cache",
                              "Error while attempting to open aerosol height file " + m_aerosol_height_filename_2);
            throw c;
        }



        // Retrieve image dimensions
        int nb_lat = p_geopotential_file->getDim("y").getSize();
        int nb_lon = p_geopotential_file->getDim("x").getSize();

        // Resize cache arrays
        m_geopotential_cache.resize(nb_lat, nb_lon);
        m_elevation_cache.resize(nb_lat, nb_lon);
        m_aerosol_height_cache_1.resize(nb_lat, nb_lon);
        m_aerosol_height_cache_2.resize(nb_lat, nb_lon);
        m_latitude_cache.resize(nb_lat, nb_lon);
        m_longitude_cache.resize(nb_lat, nb_lon);


        // Load data into cache
        logger->write_log(LogLevel::_DEBUG_, "ModelParamsReader", "update_cache",
                          "Loading cache from data files");

        p_geopotential_file->getVar("z").getVar(m_geopotential_cache.data());
        p_elevation_file->getVar("elevation").getVar(m_elevation_cache.data());
        p_aerosol_height_file_1->getVar("aerosol_height").getVar(m_aerosol_height_cache_1.data());
        p_aerosol_height_file_2->getVar("aerosol_height").getVar(m_aerosol_height_cache_2.data());

        if (!ecmwf_constant){
            p_geopotential_file->getVar("latitude").getVar(m_latitude_cache.data());
            p_geopotential_file->getVar("longitude").getVar(m_longitude_cache.data());
            m_longitude_cache = m_longitude_cache+180;
        }


        // Signal that cache was successfuly updated
        m_cache_ok = true;
        logger->write_log(LogLevel::_DEBUG_, "ModelParamsReader", "update_cache",
                          "Static data cache was successfully loaded");
    }
    else {
        logger->write_log(LogLevel::_DEBUG_, "ModelParamsReader", "update_cache",
                          "Nothing to do, static data cache already loaded");
    }
}

float ModelParamsReader::distance_geo(const float lat1, const float lon1, const float lat2, const float  lon2) {
float p = 0.017453292519943295;    // Math.PI / 180
float a = 0.5 - cos ((lat2 - lat1) * p)/2 + cos (lat1 * p) * cos (lat2 * p) *  (1 - cos ((lon2 - lon1) * p))/2;

return 12742 * asin (sqrt (a)); // 2 * R; R = 6371 km
}
