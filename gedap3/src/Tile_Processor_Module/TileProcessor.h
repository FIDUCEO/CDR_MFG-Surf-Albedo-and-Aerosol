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
#ifndef TILE_MANAGER_MODULE_TileProcessor_H_
#define TILE_MANAGER_MODULE_TileProcessor_H_
#include <random>
#include "../Common/Global_header.h"
#include "../Common/Logger.h"
#include "../Tiles/TileNetCDFManager.h"
#include "../Tile_Maker_Module/TileMaker.h"
#include "../Tile_Processor_Module/PriorManager.h"
#ifdef __cplusplus
extern "C" {
#endif
int c_cisar_cold_init(const char *lut_dir, const char *lut_id, const char *aeromod_filename);
void c_get_aerosol_coef(double *da_ext055, double *da_ext, double *da_sca, char ca_aero_model_name[][64], char *ca_aero_model_type);
void c_get_aerosol_index(const char *aerosol_class_name, int *position);
int c_cisar_inverse(Pixel_Input v_pixel_input[], Pixel_Output v_pixel_output[],
                    bool *b_temporal_smoothness, bool *b_spectral_constraint, bool *b_surface_spectral_constraint,
                    int *nb_pixels);
#ifdef __cplusplus
}
#endif
class CisarError : public std::runtime_error {
public:
 CisarError(const char *msg, int error_code)
   : std::runtime_error(std::string(msg).append(" [") + std::to_string(error_code) + "]") {}
};
class TileProcessor {
    std::string cisar_lut_directory;
    std::string cisar_lut_ID;
    std::string cisar_aeromod_filename;
    const static int min_nb_angles = 6;
    const static int min_nb_angles_per_day = 6;
    float cos_sun_glint_threshold = cos(20.0 * COEFF_DEG2RAD);
    float cos_sun_glint_affected_threshold = cos(30.0 * COEFF_DEG2RAD);
    bool b_temporal_smoothness;
    bool b_spectral_constraint;
    bool b_surface_spectral_constraint;
    int nb_pixels = 0;
    int nb_lines = 0;
    int nb_columns = 0;
    int nb_bands = 0;
    int nb_aerosols_classes = 0;
    int max_nb_iterations = 10;
    int nb_days_acc_period;
    std::vector<std::string> v_bands;
    std::vector<Aerosol_Class> v_aerosol_classes;
    int aerosol_mask[MAX_NB_AEROSOLS_CLASSES] = {0};
    double pert_first_guess = 0.4;
    double lower_FG_AOT = 1.e-14;
    double upper_FG_AOT = 0.15;
    double dust_first_guess;
    double log_a0 = 0.04878883;
    double log_a1 = 0.02386769;
    double log_a2 = -0.00116991;
    double sigmoid_width = 0.25;
    double sigmoid_max = 0.7;
    double sigmoid_center = 7.0;
    double sigmoid_offset = 0.05;
    double sigma_spectral_constraint = 1.0;
    std::shared_ptr<TileNetCDFManager> tile_netCDF_manager;
    std::shared_ptr<PriorManager> prior_manager;
    StaticData static_data;
    std::vector<InputData> v_input_data;
    std::vector<Error> v_error;
    std::vector<int> v_nb_angles;
    PriorData prior_data;
    std::vector<Pixel_Input> v_pixel_input;
    std::vector<std::vector<bool>> v_dust;
    std::vector<std::vector<int>> v_selected_angles;
    std::vector<Pixel_Output> v_pixel_output;
 TileDimensions tile_dims;
protected:
    std::shared_ptr<AbstractLogger> logger;
public:
    TileProcessor(const std::string cisar_lut_directory, const std::string cisar_lut_ID,
      const std::string cisar_aeromod_filename,
      const int max_nb_iterations, const int nb_days_acc_period, const std::vector<std::string> &v_bands,
                  const std::vector<std::string> &v_st_aerosol_classes,
      const bool b_temporal_smoothness, const bool b_spectral_constraint,
      const bool b_surface_spectral_constraint, std::shared_ptr<AbstractLogger> logger);
    virtual ~TileProcessor();
 void set_tile_dims(const TileDimensions_Tuple& tuple_tile_dims)
 {
  tile_dims = TileDimensions(std::get<0>(tuple_tile_dims), std::get<1>(tuple_tile_dims), std::get<2>(tuple_tile_dims),
      std::get<3>(tuple_tile_dims), std::get<4>(tuple_tile_dims));
  nb_lines = tile_dims.nb_lines;
  nb_columns = tile_dims.nb_columns;
  nb_pixels = tile_dims.nb_pixels;
 }
 void prepare_to_process_new_tile(const TileDimensions_Tuple& tuple_tile_dims,
   const std::map<std::string, std::vector<double>>& map_prior_rpv, const double sigma_rpv, const std::vector<double>& v_prior_AOT, const double sigma_AOT,
   const double coeff_fine_coarse_mode, double dust_prior, double dust_sigma_prior, const double dust_first_guess,
   const double pert_first_guess, const double lower_FG_AOT, const double upper_FG_AOT, const std::string aot_prior_directory)
 {
  set_tile_dims(tuple_tile_dims);
  this->pert_first_guess = pert_first_guess;
  this->lower_FG_AOT = lower_FG_AOT;
  this->upper_FG_AOT = upper_FG_AOT;
  prior_manager.reset();
  prior_manager = std::make_shared<PriorManager>(v_bands, v_aerosol_classes, nb_days_acc_period, map_prior_rpv,
                sigma_rpv, v_prior_AOT, sigma_AOT, coeff_fine_coarse_mode,
                dust_prior, dust_sigma_prior, aot_prior_directory);
 }
 void set_sigma_spectral_constraint(const double sigma)
 {
  this->sigma_spectral_constraint = sigma;
 }
 void set_params_sigmoid(const double sigmoid_width, const double sigmoid_max, const double sigmoid_center, const double sigmoid_offset)
 {
  this->sigmoid_width = sigmoid_width;
  this->sigmoid_max = sigmoid_max;
  this->sigmoid_center = sigmoid_center;
  this->sigmoid_offset = sigmoid_offset;
 }
 void set_params_log(const double log_a0, const double log_a1, const double log_a2)
 {
  this->log_a0 = log_a0;
  this->log_a1 = log_a1;
  this->log_a2 = log_a2;
 }
    void load_static_data_tile(const std::string file_static_data_tile);
 void load_input_tiles(std::vector<std::string>& v_file_input_tiles);
 void load_prior_tile(const std::string file_prior_tile);
 void create_default_prior_tile(const std::string file_prior_tile, const double start_date,
           const TileDimensions_Tuple &tuple_tile_dims,
           const std::map<std::string, std::vector<double>> &map_prior_rpv,
           const std::vector<double> &v_prior_AOT, const double sigma_rpv,
           const double sigma_AOT, const double coeff_fine_coarse_mode,
           const std::string aot_prior_directory, const std::vector<double> &latitude,
           const std::vector<double> &longitude, const bool cold_start);
    void inversion(const int period_index, const double start_acc_period_date, const double end_acc_period_date,
                   const std::string file_static_data_tile, std::vector<std::string>& v_file_input_tiles, std::vector<std::string>& satellite_file_images,
                   std::vector<std::string>& list_auxiliary_data,
                   const std::string file_prior_tile, const std::string file_next_prior_tile,
                   const std::string file_solution_tile, const bool b_training_period, const int line_stride=1, const int column_stride=1);
 std::vector<Pixel_Input> get_observations_SEVIRI(const int nb_lines, const int nb_columns, std::vector<InputData>& v_tiles, const double start_acc_period_date)
 {
  this->v_input_data = v_tiles;
        v_pixel_input.clear();
        int i_pixel = 0;
  for(int i_line = 0; i_line < nb_lines; i_line++)
  {
   for(int i_column = 0; i_column < nb_columns; i_column++)
   {
    Pixel_Input px_input = Pixel_Input(nb_bands, nb_aerosols_classes, MAX_NB_ANGLES, MAX_NB_ATM);
    std::vector<bool> v_dust_pixel(px_input.nb_obs, false);
                std::vector<int> v_selected_angles_pixel(px_input.nb_angles, -1);
                get_observations(i_line, i_column, px_input, start_acc_period_date, v_dust_pixel, v_selected_angles_pixel);
    v_pixel_input.push_back(px_input);
    i_pixel++;
   }
  }
  return v_pixel_input;
 }
private:
    void prepare_input_for_inversion(const int period_index, const double start_acc_period_date, const int line_stride,
                                     const int column_stride);
    void get_static_data(const int i_line, const int i_column, Pixel_Input &px_input);
    void
    get_observations(const int i_line, const int i_column, Pixel_Input &px_input, const double start_acc_period_date,
                     std::vector<bool> &v_dust_pixel, std::vector<int> &v_selected_angles_pixel);
    void get_prior(const int i_pixel, const int i_line, const int i_column, Pixel_Input &px_input,
                   std::vector<bool> &v_dust_pixel);
    void get_smoothness_sigma(Pixel_Input &px_input, std::vector<bool> &v_dust_pixel);
    void get_spectral_constraint_sigma(Pixel_Input &px_input);
    void get_first_guess(Pixel_Input &px_input, const int period_index, std::vector<bool> &v_dust_pixel);
};
#endif
