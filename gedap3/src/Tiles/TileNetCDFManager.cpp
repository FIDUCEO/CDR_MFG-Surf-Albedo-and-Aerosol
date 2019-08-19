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
#include "TileNetCDFManager.h"
#include <netcdf>
using std::string;
using std::runtime_error;
using std::vector;
using std::to_string;
using std::cout;
using std::endl;
using std::exception;
using netCDF::ncByte;
using netCDF::ncInt;
using netCDF::ncFloat;
using netCDF::NcDim;
using netCDF::NcGroup;
using netCDF::NcFile;
using netCDF::exceptions::NcException;
using Eigen::ArrayXf_RowMaj;
using Eigen::ArrayXXs_RowMaj;
using Eigen::ArrayXXi_RowMaj;
using Eigen::ArrayXXf_RowMaj;
using Eigen::MatrixXd_RowMaj;
void TileNetCDFManager::load_static_data(const string file, StaticData& static_data)
{
 try{
        std::unique_ptr<NcFile> tileFile1(new NcFile(file, NcFile::read));
  NcDim dimLine = tileFile1->getDim("line");
  NcDim dimColumn = tileFile1->getDim("column");
  int nb_lines = dimLine.getSize();
  int nb_columns = dimColumn.getSize();
  static_data = StaticData(nb_lines, nb_columns);
  NcGroup groupStaticData = tileFile1->getGroup("Static_Data");
  groupStaticData.getVar("Longitude").getVar(static_data.longitude.data());
  groupStaticData.getVar("Latitude").getVar(static_data.latitude.data());
  groupStaticData.getVar("LandCover").getVar(static_data.landCover.data());
  groupStaticData.getVar("Elevation").getVar(static_data.elevation.data());
    }
 catch (NcException& c){
        string error_message = "Cannot load static data tile with path " + file;
  logger->write_log(LogLevel::_ERROR_, "TileNetCDFManager", "load_static_data", error_message);
  throw runtime_error(error_message);
 }
}
void TileNetCDFManager::load_input_data(const string file, InputData& input_data)
{
try{
    std::unique_ptr<NcFile> tileFile(new NcFile(file, NcFile::read));
    string st_under_illuminated = "";
    tileFile->getAtt("is_Under_illuminated").getValues(st_under_illuminated);
    bool isUnderIlluminated = (st_under_illuminated == "TRUE");
    NcDim dimLine = tileFile->getDim("line");
    NcDim dimColumn = tileFile->getDim("column");
    NcDim dimBand = tileFile->getDim("band");
    NcDim dimDetector = tileFile->getDim("detector");
    int nb_lines = dimLine.getSize();
    int nb_columns = dimColumn.getSize();
    int nb_pixels = nb_lines * nb_columns;
    int nb_bands = dimBand.getSize();
    vector<string> v_bands;
    for(auto& group : tileFile->getGroup("Observations").getGroups()) v_bands.push_back( group.first );
    input_data = InputData(nb_lines, nb_columns, v_bands);
    NcGroup groupObservations = tileFile->getGroup("Observations");
    int i_band = 0;
    for(auto& group : groupObservations.getGroups()){
        NcGroup groupBand = group.second;
        BandObservations& band_obs_data = input_data.v_band_obs[i_band++];
        groupBand.getAtt("Channel_Wavelength").getValues(&band_obs_data.channelWaveLength);
        groupBand.getVar("TOA_BRF").getVar(band_obs_data.TOABRF.data());
        groupBand.getVar("Radiometric_Noise").getVar(band_obs_data.radiometricNoise.data());
        groupBand.getVar("Acquisition_Time").getVar(band_obs_data.acquisitionTime.data());
    }
    groupObservations.getVar("Sun_Zenith_Angle").getVar(input_data.SZA.data());
    groupObservations.getVar("Sun_Azimuth_Angle").getVar(input_data.SAA.data());
    groupObservations.getVar("Viewing_Zenith_Angle").getVar(input_data.VZA.data());
    groupObservations.getVar("Viewing_Azimuth_Angle").getVar(input_data.VAA.data());
    NcGroup groupModelParams = tileFile->getGroup("ModelParams");
    groupModelParams.getVar("TCO3").getVar(input_data.TCO3.data());
    groupModelParams.getVar("TCWV").getVar(input_data.TCWV.data());
    groupModelParams.getVar("Surface_Pressure").getVar(input_data.surfPress.data());
    groupModelParams.getVar("Wind_Speed").getVar(input_data.windSpeed.data());
    groupModelParams.getVar("Wind_Direction").getVar(input_data.windDir.data());
    groupModelParams.getVar("Aerosol_Elevation").getVar(input_data.aerosolHeight.data());
    NcGroup groupMask = tileFile->getGroup("Mask");
    if( !groupMask.isNull() ){
        groupMask.getVar("Cloud_Mask").getVar(input_data.mask.data());
        groupMask.getVar("Data_Quality_Bitmask").getVar(input_data.data_quality_bitmask.data());
    }
    }
    catch (NcException& c) {
        string error_message = "Cannot load the Input Tile with path " + file + "Image skipped.";
        logger->write_log(LogLevel::_ERROR_, "TileNetCDFManager", "load_input_data", error_message);
    }
    catch (...){
        logger->write_log(LogLevel::_ERROR_, "TileNetCDFManager", "load_input_data",
                          "Error while opening " + file + "NetCDF file. Image skipped.");
    }
}
void TileNetCDFManager::load_prior_data(const string file, PriorData& prior_data, const vector<Aerosol_Class>& v_aerosol_classes)
{
    try{
        logger->write_log(LogLevel::_DEBUG_, "TileNetCDFManager", "load_prior_data",
        "Opening NetCDF Prior Tile " + file);
        std::unique_ptr<NcFile> tileFile(new NcFile(file, NcFile::read));
  NcDim dimLine = tileFile->getDim("line");
  NcDim dimColumn = tileFile->getDim("column");
  NcDim dimBand = tileFile->getDim("band");
  NcDim dimSrfParams = tileFile->getDim("srf_param");
  NcDim dimMeanVector = tileFile->getDim("mean_vector");
  int nb_lines = dimLine.getSize();
  int nb_columns = dimColumn.getSize();
  int nb_pixels = nb_lines * nb_columns;
  int max_nb_obs = MAX_NB_OBS;
  int nb_bands = dimBand.getSize();
  vector<string> v_bands;
  for(auto& group : tileFile->getGroups()) v_bands.push_back( group.first );
  prior_data = PriorData(nb_lines, nb_columns, v_bands, v_aerosol_classes);
  tileFile->getVar("last_update_time").getVar( prior_data.last_update_time.data() );
  int band = 0;
  for(auto& group : tileFile->getGroups())
  {
   NcGroup groupBand = group.second;
   groupBand.getVar("rho_0").getVar( &prior_data.rho_0(band,0) );
   groupBand.getVar("k").getVar( &prior_data.k(band,0) );
   groupBand.getVar("theta").getVar( &prior_data.theta(band,0) );
   groupBand.getVar("rho_c").getVar( &prior_data.rho_c(band,0) );
   groupBand.getVar("sigma_rho_0").getVar( &prior_data.sigma_rho_0(band,0) );
   groupBand.getVar("sigma_k").getVar( &prior_data.sigma_k(band,0) );
   groupBand.getVar("sigma_theta").getVar( &prior_data.sigma_theta(band,0) );
   groupBand.getVar("sigma_rho_c").getVar( &prior_data.sigma_rho_c(band,0) );
   groupBand.getVar("aot_fine_abs").getVar( &prior_data.aot_fine_abs(band,0) );
   groupBand.getVar("aot_fine_non_abs").getVar( &prior_data.aot_fine_non_abs(band,0) );
   groupBand.getVar("aot_coarse").getVar( &prior_data.aot_coarse(band,0) );
   groupBand.getVar("sigma_aot").getVar( &prior_data.sigma_aot(band,0) );
   band++;
  }
  tileFile->getAtt("step_mean").getValues( &prior_data.step_mean );
  tileFile->getAtt("nb_days_since_last_reset").getValues( &prior_data.nb_days_since_last_reset );
  tileFile->getVar("nb_data").getVar( prior_data.nb_data.data() );
  tileFile->getVar("sums").getVar( &prior_data.v_sums[0] );
  tileFile->getVar("min_rpv").getVar( &prior_data.v_min_rpv[0] );
  tileFile->getVar("max_rpv").getVar( &prior_data.v_max_rpv[0] );
 }
 catch (NcException& c){
  string error_message = "Cannot load prior tile with path " + file;
  logger->write_log(LogLevel::_ERROR_, "TileNetCDFManager", "load_prior_data", error_message);
  throw c;
 }
 catch (...) {
  string error_message = "Cannot load prior tile with path " + file;
  logger->write_log(LogLevel::_ERROR_, "TileNetCDFManager", "load_prior_data", error_message);
  throw runtime_error(error_message);
 }
}
void TileNetCDFManager::save_static_data(const string file, StaticData &static_data) {
 try {
        NcFile tileFile(file, NcFile::replace, NcFile::nc4);
  NcDim dimLine = tileFile.addDim("line", static_data.nb_lines);
  NcDim dimColumn = tileFile.addDim("column", static_data.nb_columns);
  vector<NcDim> dims;
  dims.push_back(dimLine);
  dims.push_back(dimColumn);
  NcGroup groupStaticData = tileFile.addGroup("Static_Data");
        groupStaticData.addVar("Longitude", ncFloat, dims).putVar(static_data.longitude.data());
        groupStaticData.addVar("Latitude", ncFloat, dims).putVar(static_data.latitude.data());
  groupStaticData.addVar("LandCover", ncByte, dims).putVar(static_data.landCover.data());
  groupStaticData.addVar("Elevation", ncFloat, dims).putVar(static_data.elevation.data());
  logger->write_log(LogLevel::_INFO_, "TileNetCDFManager", "save_static_data",
        "Static data tile written to " + file);
 }
 catch (NcException &c) {
        string error_message = "Error while attempting to save static tile to " + file;
        logger->write_log(LogLevel::_ERROR_, "TileNetCDFManager", "save_static_data", error_message);
  throw runtime_error(error_message);
 }
}
void TileNetCDFManager::save_input_data(const string file, InputData &input_data) {
    bool isUnderIlluminated = input_data.checkIfUnderIlluminated();
    bool success = false;
    int retry = 0;
    while (!success && retry < 4) {
        try {
            std::unique_ptr<NcFile> tileFile(new NcFile(file, NcFile::replace, NcFile::nc4));
            tileFile->putAtt("is_Under_illuminated", (isUnderIlluminated) ? "TRUE" : "FALSE");
            NcDim dimLine = tileFile->addDim("line", input_data.nb_lines);
            NcDim dimColumn = tileFile->addDim("column", input_data.nb_columns);
            NcDim dimBand = tileFile->addDim("band", input_data.nb_bands);
            vector<NcDim> dims;
            dims.push_back(dimLine);
            dims.push_back(dimColumn);
            NcGroup groupObservations = tileFile->addGroup("Observations");
            int i_band = 0;
            for (const BandObservations &band_obs_data : input_data.v_band_obs) {
                NcGroup groupBand = groupObservations.addGroup(input_data.v_bands[i_band++]);
                groupBand.putAtt("Channel_Wavelength", ncFloat, band_obs_data.channelWaveLength);
                groupBand.addVar("TOA_BRF", ncFloat, dims).putVar(band_obs_data.TOABRF.data());
                groupBand.addVar("Radiometric_Noise", ncFloat, dims).putVar(band_obs_data.radiometricNoise.data());
                groupBand.addVar("Acquisition_Time", ncFloat, dims).putVar(band_obs_data.acquisitionTime.data());
            }
            groupObservations.addVar("Sun_Zenith_Angle", ncFloat, dims).putVar(input_data.SZA.data());
            groupObservations.addVar("Sun_Azimuth_Angle", ncFloat, dims).putVar(input_data.SAA.data());
            groupObservations.addVar("Viewing_Zenith_Angle", ncFloat, dims).putVar(input_data.VZA.data());
            groupObservations.addVar("Viewing_Azimuth_Angle", ncFloat, dims).putVar(input_data.VAA.data());
            NcGroup groupModelParams = tileFile->addGroup("ModelParams");
            groupModelParams.addVar("TCO3", ncFloat, dims).putVar(input_data.TCO3.data());
            groupModelParams.addVar("TCWV", ncFloat, dims).putVar(input_data.TCWV.data());
            groupModelParams.addVar("Surface_Pressure", ncFloat, dims).putVar(input_data.surfPress.data());
            groupModelParams.addVar("Wind_Speed", ncFloat, dims).putVar(input_data.windSpeed.data());
            groupModelParams.addVar("Wind_Direction", ncFloat, dims).putVar(input_data.windDir.data());
            groupModelParams.addVar("Aerosol_Elevation", ncFloat, dims).putVar(input_data.aerosolHeight.data());
            NcGroup groupMask = tileFile->addGroup("Mask");
            groupMask.addVar("Cloud_Mask", ncByte, dims).putVar(input_data.mask.data());
            groupMask.addVar("Data_Quality_Bitmask", ncByte, dims).putVar(input_data.data_quality_bitmask.data());
            logger->write_log(LogLevel::_INFO_, "TileNetCDFManager", "save_input_data",
                              "Input data tile written to " + file);
            success = true;
        }
        catch (NcException &c) {
            retry++;
        }
    }
    if (!success) {
        string error_message = "Error while attempting to save input tile to " + file;
        logger->write_log(LogLevel::_ERROR_, "TileNetCDFManager", "save_input_data", error_message);
        throw runtime_error(error_message);
    }
}
void TileNetCDFManager::save_prior_data(const string file, PriorData &prior_data) {
    try {
        std::unique_ptr<NcFile> tileFile(new NcFile(file, NcFile::replace, NcFile::nc4));
        NcDim dimLine = tileFile->addDim("line", prior_data.nb_lines);
        NcDim dimColumn = tileFile->addDim("column", prior_data.nb_columns);
        NcDim dimBand = tileFile->addDim("band", prior_data.nb_bands);
        NcDim dimAerCl = tileFile->addDim("aerosol_class", prior_data.nb_aerosol_classes);
        NcDim dimSrfParams = tileFile->addDim("srf_param", NB_SURFACE_PARAMS);
        NcDim dimMeanVector = tileFile->addDim("mean_vector", prior_data.nb_mean_vectors);
        vector<NcDim> dims_Time;
        dims_Time.push_back(dimLine);
        dims_Time.push_back(dimColumn);
        vector<NcDim> dims_Srf = dims_Time;
        vector<NcDim> dims_nb_data = {dimMeanVector, dimLine, dimColumn};
        vector<NcDim> dims_Mean{dimMeanVector, dimSrfParams, dimBand, dimLine, dimColumn};
        string list_aerosol_class_as_string = "";
        for (unsigned int i_aer = 0; i_aer < prior_data.nb_aerosol_classes; i_aer++) {
            if (i_aer > 0) list_aerosol_class_as_string += ", ";
            const Aerosol_Class &aer_class = prior_data.v_aerosol_classes[i_aer];
            string type = (aer_class.type == AEROSOL_TYPE::COARSE_MODE) ? "Coarse" : "Fine";
            list_aerosol_class_as_string += aer_class.name + " (" + type + ")";
        }
        tileFile->putAtt("list_aerosol_class", list_aerosol_class_as_string);
        string list_bands_as_string = "";
        for (unsigned int i_band = 0; i_band < prior_data.nb_bands; i_band++) {
            if (i_band > 0) list_bands_as_string += ", ";
            list_bands_as_string += prior_data.v_bands[i_band];
        }
        tileFile->putAtt("list_bands", list_bands_as_string);
        tileFile->addVar("last_update_time", ncFloat, dims_Time).putVar(prior_data.last_update_time.data());
        for (int band = 0; band < prior_data.nb_bands; band++) {
            NcGroup groupBand = tileFile->addGroup(prior_data.v_bands[band]);
            groupBand.addVar("rho_0", ncFloat, dims_Srf).putVar(&prior_data.rho_0(band, 0));
            groupBand.addVar("k", ncFloat, dims_Srf).putVar(&prior_data.k(band, 0));
            groupBand.addVar("theta", ncFloat, dims_Srf).putVar(&prior_data.theta(band, 0));
            groupBand.addVar("rho_c", ncFloat, dims_Srf).putVar(&prior_data.rho_c(band, 0));
            groupBand.addVar("sigma_rho_0", ncFloat, dims_Srf).putVar(&prior_data.sigma_rho_0(band, 0));
            groupBand.addVar("sigma_k", ncFloat, dims_Srf).putVar(&prior_data.sigma_k(band, 0));
            groupBand.addVar("sigma_theta", ncFloat, dims_Srf).putVar(&prior_data.sigma_theta(band, 0));
            groupBand.addVar("sigma_rho_c", ncFloat, dims_Srf).putVar(&prior_data.sigma_rho_c(band, 0));
            groupBand.addVar("aot_coarse", ncFloat, dims_Srf).putVar(&prior_data.aot_coarse(band, 0));
            groupBand.addVar("aot_fine_non_abs", ncFloat, dims_Srf).putVar(&prior_data.aot_fine_non_abs(band, 0));
            groupBand.addVar("aot_fine_abs", ncFloat, dims_Srf).putVar(&prior_data.aot_fine_abs(band, 0));
            groupBand.addVar("sigma_aot", ncFloat, dims_Srf).putVar(&prior_data.sigma_aot(band, 0));
        }
        tileFile->putAtt("step_mean", ncInt, prior_data.step_mean);
        tileFile->putAtt("nb_days_since_last_reset", ncInt, prior_data.nb_days_since_last_reset);
        tileFile->addVar("nb_data", ncInt, dims_nb_data).putVar(prior_data.nb_data.data());
        tileFile->addVar("sums", ncFloat, dims_Mean).putVar(&prior_data.v_sums[0]);
        tileFile->addVar("min_rpv", ncFloat, dims_Mean).putVar(&prior_data.v_min_rpv[0]);
        tileFile->addVar("max_rpv", ncFloat, dims_Mean).putVar(&prior_data.v_max_rpv[0]);
        logger->write_log(LogLevel::_INFO_, "TileNetCDFManager", "save_prior_data",
                          "Prior data tile written to" + file);
    }
    catch (NcException &c) {
        string error_message = "Error while attempting to write prior tile to " + file;
        logger->write_log(LogLevel::_ERROR_, "TileNetCDFManager", "save_prior_data", error_message);
        throw runtime_error(error_message);
    }
}
void
TileNetCDFManager::save_solution_data(const string file, vector<Pixel_Output> &v_pixel_output, vector<Pixel_Input> &v_pixel_input, vector<Error> &v_error,
                                      vector<vector<bool>> &v_dust, vector<vector<int>> &v_selected_angles,
                                      const int nb_lines, const int nb_columns, const vector<string> &v_bands,
                                      const vector<Aerosol_Class> &v_aerosol_classes,
                                      const bool b_temporal_smoothness, const bool b_spectral_constraint,
                                      const std::vector<std::string> &satellite_images_filenames,
                                      const std::vector<std::string> &list_auxiliary_data) {
    try {
        std::unique_ptr<NcFile> tileFile(new NcFile(file, NcFile::replace, NcFile::nc4));
  int nb_pixels = nb_lines * nb_columns;
  if (nb_pixels != v_pixel_output.size() || nb_pixels != v_dust.size()) {
   string error_message =
     "nb_pixels, v_pixel_output.size() and v_dust.size() must be all equal. Actual values are: " +
     to_string(nb_pixels) + ", " + to_string(v_pixel_output.size()) + ", " + to_string(v_dust.size());
   logger->write_log(LogLevel::_ERROR_, "TileNetCDFManager", "save_solution_data", error_message);
   cout << "TileNetCDFManager::save_solution_data --" << error_message << endl;
   throw runtime_error(error_message);
  }
        int max_nb_obs = 0;
        int max_nb_angles = 0;
        for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
         Pixel_Output &px_output = v_pixel_output[i_pixel];
         Pixel_Input &px_input = v_pixel_input[i_pixel];
         px_output.nb_angles = px_input.nb_angles;
            if (px_output.nb_obs > max_nb_obs) max_nb_obs = px_output.nb_obs;
            if (px_output.nb_angles > max_nb_angles) max_nb_angles = px_output.nb_angles;
        }
        NcDim dimLine = tileFile->addDim("line", nb_lines);
        NcDim dimColumn = tileFile->addDim("column", nb_columns);
        NcDim dimObs = tileFile->addDim("max_obs", max_nb_obs);
        NcDim dimObsChecked = tileFile->addDim("nr_obs_checked", v_error.size());
        vector<NcDim> dim_line_col = {dimLine, dimColumn};
        vector<NcDim> dim_obs = {dimObs, dimLine, dimColumn};
        vector<NcDim> dim_obs_checked = {dimObsChecked, dimLine, dimColumn};
#ifdef GET_ERROR
        NcDim dimAngles = tileFile->addDim("max_angles", max_nb_angles);
        vector<NcDim> dim_angles = {dimAngles, dimLine, dimColumn};
#endif
        int nb_bands = v_bands.size();
        int nb_aer_classes = v_aerosol_classes.size();
        tileFile->putAtt("temporal_smoothness_term", (b_temporal_smoothness) ? "TRUE" : "FALSE");
        tileFile->putAtt("spectral_constraint_term", (b_spectral_constraint) ? "TRUE" : "FALSE");
        string list_satellite_images_filenames_as_string = "";
        for (unsigned int i_aer = 0; i_aer < satellite_images_filenames.size(); i_aer++) {
            list_satellite_images_filenames_as_string += ", " + satellite_images_filenames[i_aer];
        }
        tileFile->putAtt("source", list_satellite_images_filenames_as_string);
        string list_auxiliary_data_as_string = "";
        for (unsigned int i_aer = 0; i_aer < list_auxiliary_data.size(); i_aer++) {
            list_auxiliary_data_as_string += ", " + list_auxiliary_data[i_aer];
        }
        tileFile->putAtt("auxiliary_data", list_auxiliary_data_as_string);
        string list_aerosol_class_as_string = "";
        for (unsigned int i_aer = 0; i_aer < v_aerosol_classes.size(); i_aer++) {
            if (i_aer > 0) list_aerosol_class_as_string += ", ";
            const Aerosol_Class &aer_class = v_aerosol_classes[i_aer];
            string type = (aer_class.type == AEROSOL_TYPE::COARSE_MODE) ? "Coarse" : "Fine";
            list_aerosol_class_as_string += aer_class.name + " (" + type + ")";
        }
        tileFile->putAtt("list_aerosol_class", list_aerosol_class_as_string);
        string list_bands_as_string = "";
        for (unsigned int i_band = 0; i_band < v_bands.size(); i_band++) {
            if (i_band > 0) list_bands_as_string += ", ";
            list_bands_as_string += v_bands[i_band];
        }
        tileFile->putAtt("list_bands", list_bands_as_string);
        int *nb_obs = new int[nb_pixels];
        int *nb_iterations = new int[nb_pixels];
        float *cost = new float[nb_pixels];
        int *error_flag = new int[nb_pixels];
        int *nb_angles = new int[nb_pixels];
        ArrayXXf_RowMaj quality_flag(max_nb_obs,nb_pixels);
        for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
            Pixel_Output &px_output = v_pixel_output[i_pixel];
            nb_obs[i_pixel] = px_output.nb_obs;
            nb_angles[i_pixel] = px_output.nb_angles;
            nb_iterations[i_pixel] = px_output.nb_iterations;
            cost[i_pixel] = px_output.cost;
            for (int i_obs = 0; i_obs < max_nb_obs; i_obs++) {
                quality_flag(i_obs,i_pixel) = px_output.quality[i_obs];
            }
            error_flag[i_pixel] = px_output.status;
        }
        tileFile->addVar("nb_obs", ncInt, dim_line_col).putVar(nb_obs);
        tileFile->addVar("nb_angles", ncInt, dim_line_col).putVar(nb_angles);
        tileFile->addVar("nb_iterations", ncInt, dim_line_col).putVar(nb_iterations);
        tileFile->addVar("cost", ncFloat, dim_line_col).putVar(cost);
        tileFile->addVar("quality_flag", ncFloat, dim_obs).putVar(quality_flag.data());
        tileFile->addVar("error_flag", ncInt, dim_line_col).putVar(error_flag);
        delete[] nb_angles;
        delete[] nb_obs;
        delete[] nb_iterations;
        delete[] cost;
        delete[] error_flag;
        NcGroup groupQInf = tileFile->addGroup("Quality_Information");
        ArrayXXf_RowMaj qual_inf_conv(max_nb_obs,nb_pixels);
        ArrayXXf_RowMaj qual_inf_bhr_outofrange(max_nb_obs,nb_pixels);
        ArrayXXf_RowMaj qual_inf_aot_outofrange(max_nb_obs,nb_pixels);
        ArrayXXf_RowMaj qual_inf_cost(max_nb_obs,nb_pixels);
        ArrayXXf_RowMaj qual_inf_jacobian(max_nb_obs,nb_pixels);
        ArrayXXf_RowMaj qual_inf_rpvEntropy(max_nb_obs,nb_pixels);
        ArrayXXf_RowMaj qual_inf_aotEntropy(max_nb_obs,nb_pixels);
        for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
            Pixel_Output &px_output = v_pixel_output[i_pixel];
            for (int i_obs = 0; i_obs < max_nb_obs; i_obs++) {
                qual_inf_conv(i_obs,i_pixel) = px_output.quality_information[0][i_obs];
                qual_inf_bhr_outofrange(i_obs,i_pixel) = px_output.quality_information[1][i_obs];
                qual_inf_aot_outofrange(i_obs,i_pixel) = px_output.quality_information[2][i_obs];
                qual_inf_cost(i_obs,i_pixel) = px_output.quality_information[3][i_obs];
                qual_inf_jacobian(i_obs,i_pixel) = px_output.quality_information[4][i_obs];
                qual_inf_rpvEntropy(i_obs,i_pixel) = px_output.quality_information[5][i_obs];
                qual_inf_aotEntropy(i_obs,i_pixel) = px_output.quality_information[6][i_obs];
            }
        }
        groupQInf.addVar("convergence", ncFloat, dim_obs).putVar(qual_inf_conv.data());
        groupQInf.addVar("bhr_out_of_range", ncFloat, dim_obs).putVar(qual_inf_bhr_outofrange.data());
        groupQInf.addVar("aot_out_of_range", ncFloat, dim_obs).putVar(qual_inf_aot_outofrange.data());
        groupQInf.addVar("cost", ncFloat, dim_obs).putVar(qual_inf_cost.data());
        groupQInf.addVar("jacobian", ncFloat, dim_obs).putVar(qual_inf_jacobian.data());
        groupQInf.addVar("rpv_entropy", ncFloat, dim_obs).putVar(qual_inf_rpvEntropy.data());
        groupQInf.addVar("aot_entropy", ncFloat, dim_obs).putVar(qual_inf_aotEntropy.data());
        {
            ArrayXXf_RowMaj acq_time(max_nb_obs,nb_pixels);
            for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
                Pixel_Output &px_output = v_pixel_output[i_pixel];
                for (int i_obs = 0; i_obs < max_nb_obs; i_obs++) {
                    acq_time(i_obs,i_pixel) = px_output.time[i_obs];
                }
            }
            tileFile->addVar("acquisition_time", ncFloat, dim_obs).putVar(acq_time.data());
        }
        {
            ArrayXXs_RowMaj dust(max_nb_obs,nb_pixels);
            for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
                for (int i_obs = 0; i_obs < max_nb_obs; i_obs++) {
                    if (i_obs < v_dust[i_pixel].size() && v_dust[i_pixel][i_obs]) {
                        dust(i_obs,i_pixel) = 1;
                    } else {
                        dust(i_obs,i_pixel) = 0;
                    }
                }
            }
            tileFile->addVar("dust", ncInt, dim_obs).putVar(dust.data());
        }
        NcGroup groupRPV = tileFile->addGroup("RPV_Parameters");
        for (int i_band = 0; i_band < nb_bands; i_band++) {
            NcGroup groupBand = groupRPV.addGroup(v_bands[i_band]);
            ArrayXf_RowMaj rho_0(nb_pixels);
            ArrayXf_RowMaj k(nb_pixels);
            ArrayXf_RowMaj theta(nb_pixels);
            ArrayXf_RowMaj rho_c(nb_pixels);
            ArrayXf_RowMaj sigma_rho_0(nb_pixels);
            ArrayXf_RowMaj sigma_k(nb_pixels);
            ArrayXf_RowMaj sigma_theta(nb_pixels);
            ArrayXf_RowMaj sigma_rho_c(nb_pixels);
            for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
                Pixel_Output &px_output = v_pixel_output[i_pixel];
                rho_0(i_pixel) = px_output.srf_param[i_band][RPV_Params_Order::RHO_0];
                k(i_pixel) = px_output.srf_param[i_band][RPV_Params_Order::K];
                theta(i_pixel) = px_output.srf_param[i_band][RPV_Params_Order::THETA];
                rho_c(i_pixel) = px_output.srf_param[i_band][RPV_Params_Order::RHO_C];
                sigma_rho_0(i_pixel) = px_output.sigma_srf_param[i_band][RPV_Params_Order::RHO_0];
                sigma_k(i_pixel) = px_output.sigma_srf_param[i_band][RPV_Params_Order::K];
                sigma_theta(i_pixel) = px_output.sigma_srf_param[i_band][RPV_Params_Order::THETA];
                sigma_rho_c(i_pixel) = px_output.sigma_srf_param[i_band][RPV_Params_Order::RHO_C];
            }
            groupBand.addVar("rho_0", ncFloat, dim_line_col).putVar(rho_0.data());
            groupBand.addVar("k", ncFloat, dim_line_col).putVar(k.data());
            groupBand.addVar("theta", ncFloat, dim_line_col).putVar(theta.data());
            groupBand.addVar("rho_c", ncFloat, dim_line_col).putVar(rho_c.data());
            groupBand.addVar("sigma_rho_0", ncFloat, dim_line_col).putVar(sigma_rho_0.data());
            groupBand.addVar("sigma_k", ncFloat, dim_line_col).putVar(sigma_k.data());
            groupBand.addVar("sigma_theta", ncFloat, dim_line_col).putVar(sigma_theta.data());
            groupBand.addVar("sigma_rho_c", ncFloat, dim_line_col).putVar(sigma_rho_c.data());
        }
        NcGroup groupAOT = tileFile->addGroup("AOT");
        for (int i_band = 0; i_band < nb_bands; i_band++) {
            NcGroup groupBand = groupAOT.addGroup(v_bands[i_band]);
            for (int i_aer = 0; i_aer < nb_aer_classes; i_aer++) {
                NcGroup groupAerClass = groupBand.addGroup(v_aerosol_classes[i_aer].name);
                ArrayXXf_RowMaj tau(max_nb_obs,nb_pixels);
                ArrayXXf_RowMaj sigma_tau(max_nb_obs,nb_pixels);
                for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
                    Pixel_Output &px_output = v_pixel_output[i_pixel];
                    for (int i_obs = 0; i_obs < max_nb_obs; i_obs++) {
                        tau(i_obs,i_pixel) = px_output.tau[i_band][i_aer][i_obs];
                        sigma_tau(i_obs,i_pixel) = px_output.sigma_tau[i_band][i_aer][i_obs];
                    }
                }
                groupAerClass.addVar("tau", ncFloat, dim_obs).putVar(tau.data());
                groupAerClass.addVar("sigma_tau", ncFloat, dim_obs).putVar(sigma_tau.data());
            }
        }
        NcGroup groupTotalAOT = tileFile->addGroup("Total_AOT");
        for (int i_band = 0; i_band < nb_bands; i_band++) {
            NcGroup groupBand = groupTotalAOT.addGroup(v_bands[i_band]);
            ArrayXXf_RowMaj total_tau(max_nb_obs,nb_pixels);
            ArrayXXf_RowMaj sigma_total_tau(max_nb_obs,nb_pixels);
            for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
                Pixel_Output &px_output = v_pixel_output[i_pixel];
                for (int i_obs = 0; i_obs < max_nb_obs; i_obs++) {
                    total_tau(i_obs,i_pixel) = px_output.total_tau[i_band][i_obs];
                    sigma_total_tau(i_obs,i_pixel) = px_output.sigma_total_tau[i_band][i_obs];
                }
            }
            groupBand.addVar("total_tau", ncFloat, dim_obs).putVar(total_tau.data());
            groupBand.addVar("sigma_total_tau", ncFloat, dim_obs).putVar(sigma_total_tau.data());
        }
        NcGroup groupAOT55 = tileFile->addGroup("AOT_550");
        for (int i_aer = 0; i_aer < nb_aer_classes; i_aer++) {
            NcGroup groupAerClass = groupAOT55.addGroup(v_aerosol_classes[i_aer].name);
            ArrayXXf_RowMaj tau_55(max_nb_obs,nb_pixels);
            ArrayXXf_RowMaj sigma_tau_55(max_nb_obs,nb_pixels);
            for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
                Pixel_Output &px_output = v_pixel_output[i_pixel];
                for (int i_obs = 0; i_obs < max_nb_obs; i_obs++) {
                    tau_55(i_obs,i_pixel) = px_output.tau_55[i_aer][i_obs];
                    sigma_tau_55(i_obs,i_pixel) = px_output.sigma_tau_55[i_aer][i_obs];
                }
            }
            groupAerClass.addVar("tau_55", ncFloat, dim_obs).putVar(tau_55.data());
            groupAerClass.addVar("sigma_tau_55", ncFloat, dim_obs).putVar(sigma_tau_55.data());
        }
        NcGroup groupAerClasProp = tileFile->addGroup("Aerosol_Class_Properties");
        for (int i_band = 0; i_band < nb_bands; i_band++) {
            NcGroup groupBand = groupAerClasProp.addGroup(v_bands[i_band]);
            ArrayXXf_RowMaj SSA(max_nb_obs,nb_pixels);
            ArrayXXf_RowMaj asym(max_nb_obs,nb_pixels);
            ArrayXXf_RowMaj sigma_SSA(max_nb_obs,nb_pixels);
            ArrayXXf_RowMaj sigma_asym(max_nb_obs,nb_pixels);
            for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
                Pixel_Output &px_output = v_pixel_output[i_pixel];
                for (int i_obs = 0; i_obs < max_nb_obs; i_obs++) {
                    SSA(i_obs,i_pixel) = px_output.SSA[i_band][i_obs];
                    sigma_SSA(i_obs,i_pixel) = px_output.sigma_SSA[i_band][i_obs];
                    asym(i_obs,i_pixel) = px_output.asym[i_band][i_obs];
                    sigma_asym(i_obs,i_pixel) = px_output.sigma_asym[i_band][i_obs];
                }
            }
            groupBand.addVar("single_scattering_albedo", ncFloat, dim_obs).putVar(SSA.data());
            groupBand.addVar("sigma_single_scattering_albedo", ncFloat, dim_obs).putVar(sigma_SSA.data());
            groupBand.addVar("asymmetry_factor", ncFloat, dim_obs).putVar(asym.data());
            groupBand.addVar("sigma_asymmetry_factor", ncFloat, dim_obs).putVar(sigma_asym.data());
        }
        NcGroup groupBRF = tileFile->addGroup("TOA_BRF");
        for (int i_band = 0; i_band < nb_bands; i_band++) {
            NcGroup groupBand = groupBRF.addGroup(v_bands[i_band]);
            ArrayXXf_RowMaj TOA_BRF(max_nb_obs,nb_pixels);
            for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
                Pixel_Output &px_output = v_pixel_output[i_pixel];
                for (int i_obs = 0; i_obs < max_nb_obs; i_obs++) {
                    TOA_BRF(i_obs,i_pixel) = px_output.TOA_BRF[i_band][i_obs];
                }
            }
            groupBand.addVar("TOA_BRF", ncFloat, dim_obs).putVar(TOA_BRF.data());
        }
#ifdef GET_ERROR
        NcGroup groupERRORS = tileFile->addGroup("ERRORS");
        for (int i_band = 0; i_band < nb_bands; i_band++) {
            NcGroup groupBand = groupERRORS.addGroup(v_bands[i_band]);
            ArrayXXf_RowMaj FASTRE_relative_err(max_nb_angles,nb_pixels);
            ArrayXXf_RowMaj EQMPN_relative_err(max_nb_angles,nb_pixels);
            ArrayXXf_RowMaj BRF_input(max_nb_angles,nb_pixels);
            ArrayXXf_RowMaj sigma_BRF_input(max_nb_angles,nb_pixels);
            ArrayXXf_RowMaj selected_angles_indices(max_nb_angles,nb_pixels);
            selected_angles_indices.setConstant(max_nb_angles,nb_pixels, -1);
            for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
                Pixel_Output &px_output = v_pixel_output[i_pixel];
                vector<int> selected_angles_indices_pixel = v_selected_angles[i_pixel];
                for (int i_angle = 0; i_angle < px_output.nb_angles; i_angle++) {
                    FASTRE_relative_err(i_angle,i_pixel) = px_output.fastre_relative_error[i_band][i_angle];
                    EQMPN_relative_err(i_angle,i_pixel) = px_output.EQMPN_relative_error[i_band][i_angle];
                    BRF_input(i_angle,i_pixel) = px_output.Input_TOA_BRF[i_band][i_angle];
                    sigma_BRF_input(i_angle,i_pixel) = px_output.Sigma_Input_TOA_BRF[i_band][i_angle];
                    selected_angles_indices(i_angle, i_pixel) = selected_angles_indices_pixel[i_angle];
                }
            }
            groupBand.addVar("FASTRE_relative_err", ncFloat, dim_angles).putVar(FASTRE_relative_err.data());
            groupBand.addVar("EQMPN_relative_err", ncFloat, dim_angles).putVar(EQMPN_relative_err.data());
            groupBand.addVar("BRF_input", ncFloat, dim_angles).putVar(BRF_input.data());
            groupBand.addVar("sigma_BRF_input", ncFloat, dim_angles).putVar(sigma_BRF_input.data());
            groupBand.addVar("selected_angles_indices", ncInt, dim_angles).putVar(selected_angles_indices.data());
        }
#endif
        NcGroup groupBHR = tileFile->addGroup("BHR");
        for (int i_band = 0; i_band < nb_bands; i_band++) {
            NcGroup groupBand = groupBHR.addGroup(v_bands[i_band]);
            ArrayXf_RowMaj BHR(nb_pixels);
            ArrayXf_RowMaj sigma_BHR(nb_pixels);
            for (int i_pixel = 0; i_pixel < nb_pixels; i_pixel++) {
                Pixel_Output &px_output = v_pixel_output[i_pixel];
                BHR(i_pixel) = px_output.BHR[i_band];
                sigma_BHR(i_pixel) = px_output.sigma_BHR[i_band];
            }
            groupBand.addVar("BHR", ncFloat, dim_line_col).putVar(BHR.data());
            groupBand.addVar("sigma_BHR", ncFloat, dim_line_col).putVar(sigma_BHR.data());
        }
        NcGroup observationValidity = tileFile->addGroup("ObservationValidity");
        ArrayXXf_RowMaj obsValidTime(v_error.size(),nb_pixels);
        ArrayXXi_RowMaj obsValidReason(v_error.size(),nb_pixels);
        for (int a = 0; a < v_error.size(); a++) {
            Error errorTemp = v_error[a];
            for (int b = 0; b < nb_lines; b++) {
                for (int c = 0; c < nb_columns; c++) {
                    obsValidTime(a, c + b * nb_columns) = errorTemp.time(b, c);
                    if (errorTemp.error[b][c] == PROCESSING_ERROR::UnderIlluminated) {
                        obsValidReason(a, c + b * nb_columns) = 1;
                    } else if (errorTemp.error[b][c] == PROCESSING_ERROR::OutOfView) {
                        obsValidReason(a, c + b * nb_columns) = 2;
                    } else if (errorTemp.error[b][c] == PROCESSING_ERROR::NegativeBRF) {
                        obsValidReason(a, c + b * nb_columns) = 5;
                    } else if (errorTemp.error[b][c] == PROCESSING_ERROR::CloudPresence) {
                        obsValidReason(a, c + b * nb_columns) = 6;
                    } else if (errorTemp.error[b][c] == PROCESSING_ERROR::SunGlint) {
                        obsValidReason(a, c + b * nb_columns) = 3;
                    } else if (errorTemp.error[b][c] == PROCESSING_ERROR::SunGlintAffected) {
                        obsValidReason(a, c + b * nb_columns) = 4;
                    } else if (errorTemp.error[b][c] == PROCESSING_ERROR::BadPixelQuality) {
                        obsValidReason(a, c + b * nb_columns) = 7;
                    }
                    else {
                        obsValidReason(a, c + b * nb_columns) = 0;
                    }
                }
            }
        }
        observationValidity.addVar("time", ncFloat, dim_obs_checked).putVar(obsValidTime.data());
        observationValidity.addVar("reasonForNotValid", ncInt, dim_obs_checked).putVar(obsValidReason.data());
        logger->write_log(LogLevel::_INFO_, "TileNetCDFManager", "save_solution_data",
                          "Solution data tile written to " + file);
    }
    catch (NcException &c) {
        string error_message = "Error while attempting to write solution tile to " + file;
        logger->write_log(LogLevel::_ERROR_, "TileNetCDFManager", "save_solution_data", error_message);
        throw runtime_error(error_message);
    }
    catch (const exception &e) {
        string error_message = "Unknown error while attempting to write solution tile to " + file +
                "\nCheck you have the right permission and that enough disk space is available."
                "\nError exception message: " + e.what();
        logger->write_log(LogLevel::_ERROR_, "TileNetCDFManager", "save_solution_data", error_message);
        throw runtime_error(error_message);
    }
}
