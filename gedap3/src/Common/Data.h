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
#ifndef TILE_MAKER_MODULE_DATA_H_
#define TILE_MAKER_MODULE_DATA_H_
#include "Global_header.h"
struct StaticData
{
 enum Land_Sea_Mask_Values {
  SPACE=0, LAND=1, SEA=2, MIXED=3, INLAND_WATER=4
 };
 int nb_lines;
 int nb_columns;
 Eigen::MatrixXf_RowMaj longitude;
 Eigen::MatrixXf_RowMaj latitude;
 Eigen::MatrixXuc_RowMaj landCover;
 Eigen::MatrixXf_RowMaj elevation;
 StaticData() : StaticData(0,0) { }
 StaticData(const int nb_lines, const int nb_columns) : nb_lines(nb_lines), nb_columns(nb_columns)
 {
  longitude = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
  latitude = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
  landCover = Eigen::MatrixXuc_RowMaj::Constant(nb_lines, nb_columns, SPACE);
  elevation = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
 }
 ~StaticData() { }
 void reset()
 {
  longitude.setConstant(MISSING_VALUE);
  latitude.setConstant(MISSING_VALUE);
  landCover.setConstant(SPACE);
  elevation.setConstant(MISSING_VALUE);
 }
 void print_data()
 {
  std::cout << std::showpoint << std::setprecision(7);
  std::cout << "longitude :" << std::endl;
  std::cout << longitude << std::endl;
  std::cout << "latitude :" << std::endl;
  std::cout << latitude << std::endl;
  std::cout << "landCover :" << std::endl;
  std::cout << landCover.cast<short>() << std::endl;
  std::cout << "elevation :" << std::endl;
  std::cout << elevation << std::endl;
 }
};
struct BandObservations
{
 float channelWaveLength;
 Eigen::MatrixXf_RowMaj TOABRF;
 Eigen::MatrixXf_RowMaj radiometricNoise;
 Eigen::MatrixXf_RowMaj acquisitionTime;
};
struct Observations
{
 int nb_bands;
 int nb_lines;
 int nb_columns;
 std::vector<BandObservations> v_band_obs;
 Eigen::MatrixXf_RowMaj SZA;
 Eigen::MatrixXf_RowMaj SAA;
 Eigen::MatrixXf_RowMaj VZA;
 Eigen::MatrixXf_RowMaj VAA;
 Observations() : Observations(0, 0, 0) { }
 Observations(const int nb_lines, const int nb_columns, const int nb_bands) :
  nb_bands(nb_bands), nb_lines(nb_lines), nb_columns(nb_columns)
 {
  for(int i_band = 0; i_band < nb_bands; i_band++)
  {
   BandObservations band_obs;
   band_obs.channelWaveLength = 0.0;
   band_obs.TOABRF = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);;
   band_obs.radiometricNoise = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);;
            band_obs.acquisitionTime = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);;
   v_band_obs.push_back(band_obs);
        }
  SZA = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
  SAA = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
  VZA = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
  VAA = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
 }
 ~Observations() { }
 void reset()
 {
  for(BandObservations& band_obs : v_band_obs)
  {
   band_obs.channelWaveLength = 0.0;
   band_obs.TOABRF.setConstant(MISSING_VALUE);
   band_obs.radiometricNoise.setConstant(MISSING_VALUE);
   band_obs.acquisitionTime.setConstant(MISSING_VALUE);
  }
  SZA.setConstant(MISSING_VALUE);
  SAA.setConstant(MISSING_VALUE);
  VZA.setConstant(MISSING_VALUE);
  VAA.setConstant(MISSING_VALUE);
 }
};
struct ModelParams
{
    enum Model_Params_Values { TC_O3, TC_WV, SURF_PRESS, WIND_SPEED, WIND_DIR, AER_HEIGHT, NB_MODEL_PARAMS };
 int nb_lines;
 int nb_columns;
 Eigen::MatrixXf_RowMaj surfPress;
 Eigen::MatrixXf_RowMaj TCO3;
 Eigen::MatrixXf_RowMaj TCWV;
 Eigen::MatrixXf_RowMaj windSpeed;
 Eigen::MatrixXf_RowMaj windDir;
    Eigen::MatrixXf_RowMaj aerosolHeight;
 ModelParams() : ModelParams(0,0) { }
 ModelParams(const int nb_lines, const int nb_columns) : nb_lines(nb_lines), nb_columns(nb_columns)
 {
  surfPress = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
  TCO3 = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
  TCWV = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
  windSpeed = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
  windDir = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
        aerosolHeight = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
 }
 ~ModelParams() { }
 void reset()
 {
  surfPress.setConstant(MISSING_VALUE);
  TCO3.setConstant(MISSING_VALUE);
  TCWV.setConstant(MISSING_VALUE);
  windSpeed.setConstant(MISSING_VALUE);
  windDir.setConstant(MISSING_VALUE);
        aerosolHeight.setConstant(MISSING_VALUE);
 }
};
struct Mask
{
    enum Mask_Values {CLOUDS, QUALITY, DUST, VOLCANIC_PLUME, SNOW, FIRE, NB_MASK_VALUES};
 int nb_lines;
 int nb_columns;
 Eigen::MatrixXuc_RowMaj mask;
    Eigen::MatrixXuc_RowMaj data_quality_bitmask;
 Mask() : Mask(0,0) { }
 Mask(const int nb_lines, const int nb_columns) : nb_lines(nb_lines), nb_columns(nb_columns)
 {
  mask = Eigen::MatrixXuc_RowMaj::Zero(nb_lines, nb_columns);
        data_quality_bitmask = Eigen::MatrixXuc_RowMaj::Zero(nb_lines, nb_columns);
 }
 ~Mask() { }
 static unsigned char bitset_to_uchar(std::bitset<NB_MASK_VALUES> value) { return (unsigned char) value.to_ulong(); }
 static std::bitset<NB_MASK_VALUES> uchar_to_bitset(unsigned char value) { return std::bitset<NB_MASK_VALUES>(value); }
 static bool check_if_cloud(unsigned char mask_val) { return std::bitset<NB_MASK_VALUES>(mask_val)[CLOUDS]; }
    static bool check_if_good(unsigned char mask_val) { return std::bitset<NB_MASK_VALUES>(mask_val)[QUALITY]; }
 static bool check_if_dust(unsigned char mask_val) { return std::bitset<NB_MASK_VALUES>(mask_val)[DUST]; }
 static bool check_if_volcanic_plume(unsigned char mask_val) { return std::bitset<NB_MASK_VALUES>(mask_val)[VOLCANIC_PLUME]; }
 static bool check_if_snow(unsigned char mask_val) { return std::bitset<NB_MASK_VALUES>(mask_val)[SNOW]; }
 static bool check_if_fire(unsigned char mask_val) { return std::bitset<NB_MASK_VALUES>(mask_val)[FIRE]; }
 void reset() { mask.setZero(); }
};
struct InputData : public Observations, public ModelParams, public Mask
{
 int nb_lines;
 int nb_columns;
 int nb_bands;
    std::vector<std::string> v_bands;
    InputData() : InputData(0, 0, std::vector<std::string>{}) {}
    InputData(const int nb_lines, const int nb_columns, const std::vector<std::string> &v_bands)
            : nb_lines(nb_lines), nb_columns(nb_columns), nb_bands(v_bands.size()), v_bands(v_bands),
              Observations(nb_lines, nb_columns, v_bands.size()), ModelParams(nb_lines, nb_columns),
              Mask(nb_lines, nb_columns) {}
 ~InputData() { }
 bool checkIfUnderIlluminated()
 {
  float minSZA = SZA.minCoeff();
  return (minSZA > SZA_THRESHOLD);
 }
};
struct PriorData
{
 int nb_lines;
 int nb_columns;
 int nb_pixels;
 int nb_bands;
    std::vector<std::string> v_bands;
 int nb_aerosol_classes;
    std::vector<Aerosol_Class> v_aerosol_classes;
 int nb_mean_vectors = 2;
 int min_data = 5;
 Eigen::VectorXd last_update_time;
 Eigen::MatrixXd_RowMaj rho_0;
 Eigen::MatrixXd_RowMaj k;
 Eigen::MatrixXd_RowMaj theta;
 Eigen::MatrixXd_RowMaj rho_c;
 Eigen::MatrixXd_RowMaj sigma_rho_0;
 Eigen::MatrixXd_RowMaj sigma_k;
 Eigen::MatrixXd_RowMaj sigma_theta;
    Eigen::MatrixXd_RowMaj sigma_rho_c;
 Eigen::MatrixXd_RowMaj aot_coarse;
 Eigen::MatrixXd_RowMaj aot_fine_non_abs;
 Eigen::MatrixXd_RowMaj aot_fine_abs;
 Eigen::MatrixXd_RowMaj sigma_aot;
    int dim_vector = 0;
    int nb_days_reset = 45;
 int nb_days_since_last_reset = 0;
 int step_mean = 0;
 Eigen::MatrixXi_RowMaj nb_data;
    Eigen::MatrixXi_RowMaj quality;
 std::vector<float> v_sums;
 std::vector<float> v_max_rpv;
 std::vector<float> v_min_rpv;
    PriorData() : PriorData(0, 0, std::vector<std::string> {}, std::vector<Aerosol_Class> {}) {}
    PriorData(const int nb_lines, const int nb_columns, const std::vector<std::string> &v_bands,
              const std::vector<Aerosol_Class> &v_aerosol_classes) :
            nb_lines(nb_lines), nb_columns(nb_columns), nb_pixels(nb_lines * nb_columns), nb_bands(v_bands.size()),
            v_bands(v_bands), nb_aerosol_classes(v_aerosol_classes.size()), v_aerosol_classes(v_aerosol_classes) {
        last_update_time = Eigen::VectorXd::Constant(nb_pixels, MISSING_VALUE);
        rho_0 = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_pixels, MISSING_VALUE);
        k = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_pixels, MISSING_VALUE);
        theta = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_pixels, MISSING_VALUE);
        rho_c = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_pixels, MISSING_VALUE);
        sigma_rho_0 = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_pixels, MISSING_VALUE);
        sigma_k = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_pixels, MISSING_VALUE);
        sigma_theta = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_pixels, MISSING_VALUE);
        sigma_rho_c = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_pixels, MISSING_VALUE);
        aot_coarse = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_pixels, MISSING_VALUE);
        aot_fine_non_abs = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_pixels, MISSING_VALUE);
        aot_fine_abs = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_pixels, MISSING_VALUE);
        sigma_aot = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_pixels, MISSING_VALUE);
        nb_data = Eigen::MatrixXi_RowMaj::Zero(nb_mean_vectors, nb_pixels);
        quality = Eigen::MatrixXi_RowMaj::Zero(nb_mean_vectors, nb_pixels);
        dim_vector = nb_mean_vectors * NB_SURFACE_PARAMS * nb_bands * nb_pixels;
        v_sums = std::vector<float>(dim_vector, 0.0);
        v_max_rpv = std::vector<float>(dim_vector, 0.0);
        v_min_rpv = std::vector<float>(dim_vector, 0.0);
    }
 float& get_sum(const int i_vector, const int i_rpv_param, const int i_band, const int i_pixel)
 {
  int index = get_index(i_vector, i_rpv_param, i_band, i_pixel);
        return v_sums[index];
 }
 float& get_max_rpv(const int i_vector, const int i_rpv_param, const int i_band, const int i_pixel)
 {
  int index = get_index(i_vector, i_rpv_param, i_band, i_pixel);
        return v_max_rpv[index];
 }
 float& get_min_rpv(const int i_vector, const int i_rpv_param, const int i_band, const int i_pixel)
 {
  int index = get_index(i_vector, i_rpv_param, i_band, i_pixel);
        return v_min_rpv[index];
 }
private:
   int get_index(const int i_vector, const int i_rpv_param, const int i_band, const int i_pixel)
   {
     return ((i_vector * NB_SURFACE_PARAMS + i_rpv_param) * nb_bands + i_band) * nb_pixels + i_pixel;
   }
};
#endif
