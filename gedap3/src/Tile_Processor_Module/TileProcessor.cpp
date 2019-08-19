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
#include "TileProcessor.h"
using namespace std;
TileProcessor::TileProcessor(const string cisar_lut_directory, const string cisar_lut_ID,
                             const string cisar_aeromod_filename,
                             const int max_nb_iterations, const int nb_days_acc_period, const vector<string> &v_bands,
                             const vector<string> &v_st_aerosol_classes,
                             const bool b_temporal_smoothness, const bool b_spectral_constraint,
                             const bool b_surface_spectral_constraint, shared_ptr<AbstractLogger> logger)
{
    this->cisar_lut_directory = cisar_lut_directory;
    this->cisar_lut_ID = cisar_lut_ID;
    this->cisar_aeromod_filename = cisar_aeromod_filename;
    this->logger = logger;
    this->b_temporal_smoothness = b_temporal_smoothness;
    this->b_spectral_constraint = b_spectral_constraint;
    this->b_surface_spectral_constraint = b_surface_spectral_constraint;
    this->dust_first_guess = dust_first_guess;
    if (max_nb_iterations > 0) this->max_nb_iterations = max_nb_iterations;
    this->nb_days_acc_period = nb_days_acc_period;
    this->v_bands = v_bands;
    nb_bands = v_bands.size();
    clock_t start = clock();
    int cisar_status =
            c_cisar_cold_init(cisar_lut_directory.c_str(), cisar_lut_ID.c_str(), cisar_aeromod_filename.c_str());
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "TileProcessor",
                      string("Time to call c_cisar_cold_init: ").append(
                          to_string((double) (clock() - start) / (double) CLOCKS_PER_SEC)));
    if (cisar_status > 1000) {
        logger->write_log(LogLevel::_ERROR_, "TileProcessor", "TileProcessor",
                          string("Fatal error during cold initialisation [" + to_string(cisar_status) + "]"));
        throw CisarError("Fatal error during cold initialisation", cisar_status);
    }
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "TileProcessor", "End of CISAR coldinit");
    int position=-1;
    double da_ext055[MAX_NB_AEROSOLS_CLASSES];
    double da_ext[MAX_NB_BANDS][MAX_NB_AEROSOLS_CLASSES];
    double da_sca[MAX_NB_BANDS][MAX_NB_AEROSOLS_CLASSES];
    char ca_aero_model_name[MAX_NB_AEROSOLS_CLASSES][64];
    char ca_aero_model_type[MAX_NB_AEROSOLS_CLASSES];
    c_get_aerosol_coef(&da_ext055[0], &da_ext[0][0], &da_sca[0][0], &ca_aero_model_name[0], &ca_aero_model_type[0]);
    for (string p_aer_class : v_st_aerosol_classes)
    {
        string name = p_aer_class;
        c_get_aerosol_index(name.c_str(), &position);
        if(position==-1){
            string error_message = "Aerosol Class "+name+" not found.";
            logger->write_log(LogLevel::_ERROR_, "TileProcessor", "inversion",
                              error_message);
            throw runtime_error(error_message);
        }
        AEROSOL_TYPE type = (ca_aero_model_type[position]=='C') ? AEROSOL_TYPE::COARSE_MODE : AEROSOL_TYPE::FINE_MODE;
        this->v_aerosol_classes.push_back(Aerosol_Class(name, type));
    }
    nb_aerosols_classes = this->v_aerosol_classes.size();
    position = -1;
    for (Aerosol_Class &aer_class : v_aerosol_classes)
    {
        c_get_aerosol_index(aer_class.name.c_str(), &position);
        if (position >= 0 && position < MAX_NB_AEROSOLS_CLASSES)
            aerosol_mask[position] = 1;
    }
    tile_netCDF_manager = make_shared<TileNetCDFManager>(logger);
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "constructor",
                      "An instance of TileProcessor has been created");
}
TileProcessor::~TileProcessor()
{
    v_pixel_input.clear();
    v_pixel_output.clear();
    v_input_data.clear();
    v_dust.clear();
    v_selected_angles.clear();
}
void TileProcessor::load_static_data_tile(const string file_static_data_tile)
{
    tile_netCDF_manager->load_static_data(file_static_data_tile, static_data);
    if (static_data.nb_lines != nb_lines or static_data.nb_columns != nb_columns)
    {
        string error_message = "Dimensions of current static data tile " + file_static_data_tile +
                " is not equal to the expected tile dimensions.\n";
        error_message +=
                "Current tile:  (" + to_string(static_data.nb_lines) + ", " + to_string(static_data.nb_columns) + ").\n";
        error_message += "Expected dimensions: (" + to_string(nb_lines) + ", " + to_string(nb_columns) + ")";
        logger->write_log(LogLevel::_ERROR_, "TileProcessor", "load_static_data", error_message);
        throw runtime_error(error_message);
    }
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "load_static_data",
                      "Static Tile " + file_static_data_tile + " has been loaded");
}
void TileProcessor::load_prior_tile(const string file_prior_tile)
{
    tile_netCDF_manager->load_prior_data(file_prior_tile, prior_data, v_aerosol_classes);
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "load_prior_data",
                      "Prior Tile " + file_prior_tile + " has been loaded");
}
void TileProcessor::load_input_tiles(vector<string> &v_file_input_tiles)
{
    for (string file_input_tile : v_file_input_tiles)
    {
        InputData input_data;
        tile_netCDF_manager->load_input_data(file_input_tile, input_data);
        Error error;
        error = Error(nb_lines, nb_columns);
        if (input_data.nb_lines != nb_lines or input_data.nb_columns != nb_columns)
        {
            string error_message =
                    "Dimensions of input tile " + file_input_tile + " is not equal to the expected tile dimensions.\n";
            error_message +=
                    "Current tile:  (" + to_string(input_data.nb_lines) + ", " + to_string(input_data.nb_columns) + ").\n";
            error_message += "Expected dimensions: (" + to_string(nb_lines) + ", " + to_string(nb_columns) + ")";
            logger->write_log(LogLevel::_ERROR_, "TileProcessor", "load_input_tiles", error_message);
            continue;
        }
        v_input_data.push_back(input_data);
        v_error.push_back(error);
    }
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "load_input_tiles", "All input tiles have been loaded");
}
void TileProcessor::create_default_prior_tile(const string file_prior_tile, const double start_date, const TileDimensions_Tuple& tuple_tile_dims,
                                              const map<string, vector<double>>& map_prior_rpv, const vector<double>& v_prior_AOT, const double sigma_rpv, const double sigma_AOT,
                                              const double coeff_fine_coarse_mode, const string aot_prior_directory, const vector<double>& latitude, const vector<double>& longitude, const bool cold_start)
{
    set_tile_dims(tuple_tile_dims);
    prior_manager.reset();
    prior_manager = make_shared<PriorManager>(v_bands, v_aerosol_classes, nb_days_acc_period, map_prior_rpv, sigma_rpv,
                                              v_prior_AOT, sigma_AOT, coeff_fine_coarse_mode, 0.0, 0.0,
                                              aot_prior_directory);
    if (!ifstream((file_prior_tile).c_str()).good())
    {
        v_pixel_input.clear();
        for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++)
        {
            Pixel_Input px_input(nb_bands, nb_aerosols_classes, MAX_NB_ANGLES, MAX_NB_ATM);
            v_pixel_input.push_back(px_input);
        }
        try
        {
            prior_data = prior_manager->create_default_prior(nb_lines, nb_columns, start_date, latitude, longitude,
                                                             cisar_lut_directory, cisar_lut_ID);
            tile_netCDF_manager->save_prior_data(file_prior_tile, prior_data);
            logger->write_log(LogLevel::_INFO_, "TileProcessor", "create_default_prior",
                              "Prior Tile " + file_prior_tile + " has been created with default values");
        }
        catch (runtime_error &e)
        {
            prior_manager.reset();
            string error_message = "Runtime error thrown. No Tile with default prior created";
            logger->write_log(LogLevel::_ERROR_, "TileProcessor", "create_default_prior", error_message);
            throw e;
        }
        v_pixel_input.clear();
        logger->write_log(LogLevel::_INFO_, "TileProcessor", "create_default_prior",
                          "New Prior Tile " + file_prior_tile + " created");
    }
    else
    {
        logger->write_log(LogLevel::_INFO_, "TileProcessor", "create_default_prior",
                          "Prior Tile " + file_prior_tile + " already exists. No new netCDF file created");
    }
    prior_manager.reset();
}
void TileProcessor::inversion(const int period_index, const double start_acc_period_date, const double end_acc_period_date,
                              const string file_static_data_tile, vector<string>& v_file_input_tiles, vector<string> &satellite_images_filenames,
                              vector<string> &list_auxiliary_data, const string file_prior_tile, const string file_next_prior_tile,
                              const string file_solution_tile, const bool b_training_period, const int line_stride, const int column_stride)
{
    string message = "";
    clock_t start;
    double duration_cisar_init = 0.0;
    double duration_prepare_input = 0.0;
    double duration_inversion = 0.0;
    double duration_update_prior = 0.0;
    double duration_save_output = 0.0;
    nb_days_acc_period = static_cast<int>(end_acc_period_date - start_acc_period_date) + 1;
    prior_manager->set_nb_days_acc_period(nb_days_acc_period);
    cout << showpoint << setprecision(7) << "start_acc_period_date, end_acc_period_date, nb_days_acc_period: "
         << start_acc_period_date << " " << end_acc_period_date << " " << nb_days_acc_period << endl;
    message = "Start a new inversion for period number " + to_string(period_index);
    cout << message << endl;
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "inversion", message);
    load_static_data_tile(file_static_data_tile);
    load_input_tiles(v_file_input_tiles);
    load_prior_tile(file_prior_tile);
    start = clock();
    prepare_input_for_inversion(period_index, start_acc_period_date, line_stride, column_stride);
    duration_prepare_input = (double) (clock() - start);
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "inversion", "End of prepare_input_for_inversion");
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "inversion", "Start of CISAR inversion");
    start = clock();
    int cisar_status = c_cisar_inverse(v_pixel_input.data(), v_pixel_output.data(), &b_temporal_smoothness,
                                       &b_spectral_constraint, &b_surface_spectral_constraint, &nb_pixels);
    duration_inversion = (double) (clock() - start);
    if (cisar_status > 1000) {
        CisarError e("Fatal error during inversion", cisar_status);
        logger->write_log(LogLevel::_ERROR_, "TileProcessor", "inversion",
                          e.what());
        throw e;
    }
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "inversion", "End of CISAR inversion");
    start = clock();
    tile_netCDF_manager->save_solution_data(file_solution_tile, v_pixel_output, v_pixel_input, v_error, v_dust, v_selected_angles, nb_lines, nb_columns,
                                            v_bands, v_aerosol_classes, b_temporal_smoothness, b_spectral_constraint, satellite_images_filenames, list_auxiliary_data);
    duration_save_output = (double) (clock() - start);
    v_error.clear();
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "inversion",
                      "Solution Tile " + file_solution_tile + " created");
    start = clock();
    prior_manager->update_prior(v_pixel_input, v_pixel_output, end_acc_period_date, prior_data);
    tile_netCDF_manager->save_prior_data(file_next_prior_tile, prior_data);
    duration_update_prior = (double) (clock() - start);
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "inversion",
                      "End of update_prior. Next Prior Tile " + file_next_prior_tile + " created.");
    v_pixel_input.clear();
    v_pixel_output.clear();
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "inversion",
                      string("Comp. time to prepare input vector:").append(
                          to_string(duration_prepare_input / (double) CLOCKS_PER_SEC)));
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "inversion",
                      string("Comp. time to perform inversion:").append(
                          to_string(duration_inversion / (double) CLOCKS_PER_SEC)));
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "inversion",
                      string("Comp. time to update the prior:").append(
                          to_string(duration_update_prior / (double) CLOCKS_PER_SEC)));
    logger->write_log(LogLevel::_INFO_, "TileProcessor", "inversion",
                      string("Comp. time to save the solution:").append(
                          to_string(duration_save_output / (double) CLOCKS_PER_SEC)));
}
void TileProcessor::prepare_input_for_inversion(const int period_index, const double start_acc_period_date,
                                                const int line_stride, const int column_stride)
{
    clock_t start = clock();
    v_pixel_input.clear();
    v_pixel_output.clear();
    v_dust.clear();
    v_selected_angles.clear();
    v_pixel_input.reserve(nb_pixels);
    v_pixel_output.reserve(nb_pixels);
    v_dust.reserve(nb_pixels);
    v_selected_angles.reserve(nb_pixels);
    int i_pixel = 0;
    for (int i_line = 0; i_line < nb_lines; i_line += line_stride)
    {
        for (int i_column = 0; i_column < nb_columns; i_column += column_stride)
        {
            i_pixel = i_line * nb_columns + i_column;
            vector<bool> v_dust_pixel;
            vector<int> v_selected_angles_pixel;
            Pixel_Input px_input(nb_bands, nb_aerosols_classes, MAX_NB_ANGLES, MAX_NB_ATM);
            Pixel_Output px_output(nb_bands, nb_aerosols_classes, MAX_NB_ANGLES, MAX_NB_ATM);
            px_input.max_nb_iterations = max_nb_iterations;
            get_static_data(i_line, i_column, px_input);
            get_observations(i_line, i_column, px_input, start_acc_period_date, v_dust_pixel, v_selected_angles_pixel);
            get_prior(i_pixel, i_line, i_column, px_input, v_dust_pixel);
            if (b_temporal_smoothness) get_smoothness_sigma(px_input, v_dust_pixel);
            if (b_spectral_constraint) get_spectral_constraint_sigma(px_input);
            get_first_guess(px_input, period_index, v_dust_pixel);
            v_pixel_input.push_back(px_input);
            v_pixel_output.push_back(px_output);
            v_dust.push_back(v_dust_pixel);
            v_selected_angles.push_back(v_selected_angles_pixel);
        }
    }
    v_input_data.clear();
}
void TileProcessor::get_static_data(const int i_line, const int i_column, Pixel_Input &px_input)
{
    px_input.srf_type = (int) static_data.landCover(i_line, i_column);
    if (px_input.srf_type == StaticData::MIXED) px_input.srf_type = StaticData::LAND;
    if (px_input.srf_type == StaticData::INLAND_WATER) px_input.srf_type = StaticData::SEA;
    for (int i_aer = 0; i_aer < MAX_NB_AEROSOLS_CLASSES; i_aer++) {
        px_input.aerosol_mask[i_aer] = aerosol_mask[i_aer];
    }
}
void TileProcessor::get_observations(const int i_line, const int i_column, Pixel_Input &px_input,
                                     const double start_acc_period_date, vector<bool> &v_dust_pixel, vector<int> &v_selected_angles_pixel)
{
    vector<int> v_selected_angle_indexes;
    v_selected_angle_indexes.reserve(v_input_data.size());
    vector<bool> v_dust_all_obs;
    v_dust_all_obs.reserve(v_input_data.size());
    bool b_angles_in_valid_range = false;
    bool b_brf_gt_zero = false;
    bool b_cloud_free = false;
    bool b_no_sun_glint = false;
    bool b_sun_glint_affected = false;
    bool b_quality = false;
    double mod_jul_date = 0.0;
    int start_day = -1;
    int year = -1, day = -1, dayOfTheYear = -1, hour = -1, minutes = -1;
    const double eps_time = 0.0005;
    vector<double> v_min_brf_error(24 * nb_days_acc_period, 100.0);
    vector<double> v_index_of_min_brf_error(24 * nb_days_acc_period, -1);
    vector<vector<int>> v_valid_obs_indexes;
    for (int day = 0; day < nb_days_acc_period; day++) v_valid_obs_indexes.push_back(vector<int>());
    for (int i_angle = 0; i_angle < v_input_data.size(); i_angle++)
    {
        InputData &input_data = v_input_data[i_angle];
        b_angles_in_valid_range =
                input_data.SZA(i_line, i_column) > MIN_SZA && input_data.SZA(i_line, i_column) < MAX_SZA &&
                input_data.VZA(i_line, i_column) > MIN_VZA && input_data.VZA(i_line, i_column) < MAX_VZA;
        b_brf_gt_zero = true;
        for (const BandObservations &band_obs_data : input_data.v_band_obs)
            b_brf_gt_zero &= band_obs_data.TOABRF(i_line, i_column) > EPSILON;
        b_cloud_free = !Mask::check_if_cloud(input_data.mask(i_line, i_column));
        v_dust_all_obs.push_back(Mask::check_if_dust(input_data.mask(i_line, i_column)));
        /*only use data where data_quality flag is 0 or 4 where 4 is the suspicious space cound (which is mitigated; see FCDR ATBD)*/
        b_quality = (input_data.data_quality_bitmask(i_line, i_column)==0)||(input_data.data_quality_bitmask(i_line, i_column)==4);
        if (px_input.srf_type == SEA)
        {
            double SZA = input_data.SZA(i_line, i_column) * COEFF_DEG2RAD;
            double VZA = input_data.VZA(i_line, i_column) * COEFF_DEG2RAD;
            double RAA = (input_data.VAA(i_line, i_column) - input_data.SAA(i_line, i_column)) * COEFF_DEG2RAD;
            if (RAA < 0.0) RAA += TWO_PI;
            if (RAA > PI) RAA = TWO_PI - RAA;
            double cosg = cos(SZA) * cos(VZA) - sin(SZA) * sin(VZA) * cos(RAA);
            b_no_sun_glint = cosg < cos_sun_glint_threshold;
            b_sun_glint_affected = cosg > cos_sun_glint_affected_threshold;
        }
        else
        {
            b_no_sun_glint = true;
        }
        if (b_angles_in_valid_range && b_brf_gt_zero && b_cloud_free && b_no_sun_glint && b_quality)
        {
            mod_jul_date = input_data.v_band_obs[0].acquisitionTime(i_line, i_column) + eps_time;
            MJD_to_Time(mod_jul_date, &year, &dayOfTheYear, &hour, &minutes);
            day = floor(mod_jul_date) - floor(start_acc_period_date + eps_time);
            hour = day * 24 + hour;
            double sum_brf_error = 0.0;
            for (const BandObservations &band_obs_data : input_data.v_band_obs)
                sum_brf_error += band_obs_data.radiometricNoise(i_line, i_column);
            if (sum_brf_error < v_min_brf_error[hour])
            {
                v_min_brf_error[hour] = sum_brf_error;
                v_index_of_min_brf_error[hour] = i_angle;
            }
            v_valid_obs_indexes[day].push_back(i_angle);
        }
        else
        {
        }
        if (!b_angles_in_valid_range)
        {
            v_error[i_angle].time(i_line, i_column) = input_data.v_band_obs[0].acquisitionTime(i_line, i_column);
            if (!(input_data.SZA(i_line, i_column) > MIN_SZA && input_data.SZA(i_line, i_column) < MAX_SZA))
            {
                v_error[i_angle].error[i_line][i_column] = PROCESSING_ERROR::UnderIlluminated;
            }
            else
            {
                v_error[i_angle].error[i_line][i_column] = PROCESSING_ERROR::OutOfView;
            }
        }
        else if (!b_brf_gt_zero)
        {
            v_error[i_angle].error[i_line][i_column] = PROCESSING_ERROR::NegativeBRF;
        }
        else if (!b_cloud_free)
        {
            v_error[i_angle].error[i_line][i_column] = PROCESSING_ERROR::CloudPresence;
        }
        else if (!b_no_sun_glint)
        {
            v_error[i_angle].error[i_line][i_column] = PROCESSING_ERROR::SunGlint;
        }
        else if (b_sun_glint_affected)
        {
            v_error[i_angle].error[i_line][i_column] = PROCESSING_ERROR::SunGlintAffected;
        }
        else if (b_quality)
        {
            v_error[i_angle].error[i_line][i_column] = PROCESSING_ERROR::BadPixelQuality;
        }
        else
        {
            v_error[i_angle].error[i_line][i_column] = PROCESSING_ERROR::None;
        }
    }
    for (day = 0; day < nb_days_acc_period; day++)
    {
        v_selected_angle_indexes.insert(v_selected_angle_indexes.end(), v_valid_obs_indexes[day].begin(),
                                        v_valid_obs_indexes[day].end());
    }
    if (v_selected_angle_indexes.size() < min_nb_angles || v_selected_angle_indexes.size() > MAX_NB_ANGLES)
        v_selected_angle_indexes.clear();
    v_dust_pixel.clear();
    v_dust_pixel.push_back(false);
    px_input.SrfPres[0] = 0.0;
    px_input.WV_COL[0] = 0.0;
    px_input.OZ_COL[0] = 0.0;
    px_input.layer_height[0] = 0.0;
    px_input.wind_speed = 0.0;
    px_input.wind_direction = 0.0;
    double tmp_wind_direction = 0.0;
    int i_angle = 0;
    for (int &i_selected_angle: v_selected_angle_indexes)
    {
        InputData &input_data = v_input_data[i_selected_angle];
        v_dust_pixel[0] = v_dust_pixel[0] || v_dust_all_obs[i_selected_angle];
        px_input.cloud_mask[i_angle] = input_data.mask(i_line, i_column);
        px_input.time[i_angle] = input_data.v_band_obs[0].acquisitionTime(i_line, i_column);
        px_input.SZA[i_angle] = input_data.SZA(i_line, i_column) * COEFF_DEG2RAD;
        px_input.VZA[i_angle] = input_data.VZA(i_line, i_column) * COEFF_DEG2RAD;
        px_input.RAA[i_angle] = (input_data.VAA(i_line, i_column) - input_data.SAA(i_line, i_column)) * COEFF_DEG2RAD;
        if (px_input.RAA[i_angle] < 0.0) px_input.RAA[i_angle] += TWO_PI;
        if (px_input.RAA[i_angle] > PI) px_input.RAA[i_angle] = TWO_PI - px_input.RAA[i_angle];
        px_input.SrfPres[0] += input_data.surfPress(i_line, i_column);
        px_input.WV_COL[0] += input_data.TCWV(i_line, i_column);
        px_input.OZ_COL[0] += input_data.TCO3(i_line, i_column);
        px_input.layer_height[0] += input_data.aerosolHeight(i_line, i_column);
        px_input.wind_speed += input_data.windSpeed(i_line, i_column);
        tmp_wind_direction = input_data.windDir(i_line, i_column);
        if (tmp_wind_direction > PI) tmp_wind_direction = TWO_PI - tmp_wind_direction;
        px_input.wind_direction += tmp_wind_direction;
        for (int i_band = 0; i_band < nb_bands; i_band++)
        {
            BandObservations &band_obs_data = input_data.v_band_obs[i_band];
            px_input.TOA_BRF[i_band][i_angle] = band_obs_data.TOABRF(i_line, i_column);
            px_input.sigma_TOA_BRF[i_band][i_angle] = band_obs_data.radiometricNoise(i_line, i_column);
        }
        i_angle++;
    }
    px_input.nb_angles = i_angle;
    if (px_input.nb_angles==0) px_input.nb_obs = 0;
    else px_input.nb_obs = 1;
    px_input.SrfPres[0] /= (double) px_input.nb_angles;
    px_input.WV_COL[0] /= (double) px_input.nb_angles;
    px_input.OZ_COL[0] /= (double) px_input.nb_angles;
    px_input.layer_height[0] /= (double) px_input.nb_angles;
    if (px_input.nb_angles==0){
        px_input.wind_direction = 0.0;
        px_input.wind_speed = 0.0;
    }
    else{
        px_input.wind_speed /= (double) px_input.nb_angles;
        px_input.wind_direction = 0.0;
    }
    v_selected_angles_pixel =v_selected_angle_indexes;
}
void TileProcessor::get_prior(const int i_pixel, const int i_line, const int i_column, Pixel_Input &px_input,
                              vector<bool> &v_dust_pixel)
{
    prior_manager->set_prior_to_fortran_input_structure(i_pixel, static_data.latitude(i_line, i_column),
                                                        static_data.longitude(i_line, i_column), px_input,
                                                        prior_data.last_update_time[i_pixel], prior_data, v_dust_pixel,
                                                        cisar_lut_directory, cisar_lut_ID);
}
void TileProcessor::get_smoothness_sigma(Pixel_Input &px_input, vector<bool> &v_dust_pixel)
{
    double lag = 0.0;
    double sigma = 0.0;
    for (int i_obs = 0; i_obs < px_input.nb_obs - 1; i_obs++)
    {
        lag = (px_input.time[i_obs + 1] - px_input.time[i_obs]) * 24.0;
        sigma = sigmoid_offset + (sigmoid_max / (1.0 + exp(-sigmoid_width * (lag - sigmoid_center)))) * (lag / 24.0);
        for (int i_band = 0; i_band < nb_bands; i_band++) px_input.sigma_temp_smoothness[i_band][i_obs] = sigma;
    }
}
void TileProcessor::get_spectral_constraint_sigma(Pixel_Input &px_input)
{
    px_input.sigma_spectral_constraint = sigma_spectral_constraint;
}
void TileProcessor::get_first_guess(Pixel_Input &px_input, const int period_index, vector<bool> &v_dust_pixel)
{
    double sign_srf = pow(-1.0, period_index);
    for (int i_band = 0; i_band < px_input.nb_bands; i_band++)
    {
        if (px_input.srf_type == LAND)
        {
            for (int i_rpv_param = 0; i_rpv_param < NB_SURFACE_PARAMS; i_rpv_param++)
            {
                if (px_input.sigma_prior_srf[i_band][i_rpv_param] < 0.5 * px_input.prior_srf[i_band][i_rpv_param])
                {
                    px_input.first_guess_srf[i_band][i_rpv_param] = px_input.prior_srf[i_band][i_rpv_param] +
                            sign_srf *
                            px_input.sigma_prior_srf[i_band][i_rpv_param];
                }
                else
                {
                    px_input.first_guess_srf[i_band][i_rpv_param] = px_input.prior_srf[i_band][i_rpv_param] +
                            sign_srf * 0.5 *
                            px_input.prior_srf[i_band][i_rpv_param];
                }
                check_RPV_physical_range(i_rpv_param, px_input.first_guess_srf[i_band][i_rpv_param]);
            }
        }
        else
        {
            for (int i_rpv_param = 0; i_rpv_param < NB_SURFACE_PARAMS; i_rpv_param++)
            {
                px_input.first_guess_srf[i_band][i_rpv_param] = 0.0;
            }
        }
        for (int i_aer = 0; i_aer < px_input.nb_aerosols; i_aer++)
        {
            for (int i_obs = 0; i_obs < px_input.nb_obs; i_obs++)
            {
                if (v_dust_pixel[i_obs] && v_aerosol_classes[i_aer].type == AEROSOL_TYPE::COARSE_MODE)
                {
                    px_input.first_guess_tau[i_band][i_aer][i_obs] = dust_first_guess;
                }
                else
                {
                    px_input.first_guess_tau[i_band][i_aer][i_obs] = (i_obs % 2 == 0) ? lower_FG_AOT : upper_FG_AOT;
                }
            }
        }
    }
}
