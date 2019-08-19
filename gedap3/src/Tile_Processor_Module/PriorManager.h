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
#ifndef SOURCE_CPP_SOURCE_TILE_PROCESSOR_MODULE_PRIORMANAGER_H_
#define SOURCE_CPP_SOURCE_TILE_PROCESSOR_MODULE_PRIORMANAGER_H_
#include <netcdf>
#include "../Common/Global_header.h"
#include "../Common/Logger.h"
#include "../Common/Data.h"
#include "../Tiles/FortranStructures.h"
#ifdef __cplusplus
extern "C" {
#endif
void c_get_aerosol_coef(double *da_ext055, double *da_ext, double *da_sca, char ca_aero_model_name[][64], char *ca_aero_model_type);
void c_get_aerosol_index(const char *aerosol_class_name, int *position);
#ifdef __cplusplus
}
#endif
struct PriorManager
{
 int nb_bands;
 std::vector<std::string> v_bands;
 int nb_aerosol_classes;
 std::vector<Aerosol_Class> v_aerosol_classes;
 int nb_days_acc_period;
 Eigen::MatrixXd_RowMaj prior_AOT;
 Eigen::MatrixXd_RowMaj prior_AOT_coarse;
 Eigen::MatrixXd_RowMaj prior_AOT_fine_total;
 Eigen::MatrixXd_RowMaj prior_AOT_fine_abs;
 float latitudes[180];
 float longitudes[360];
 Eigen::VectorXd aerosol_class_weights;
 double sigma_AOT;
 Eigen::MatrixXd_RowMaj default_RPV_values;
 Eigen::Vector4d default_RPV_sigma_values;
 double factor_sigma = 0.1;
 double min_threshold_sigma = 0.05;
 double max_threshold_sigma = 0.4;
 double max_daily_change = 0.02;
 double default_increase_sigma_factor_per_day = 1.05;
 double default_increase_sigma_factor = 0.0;
 double dust_prior;
 double dust_sigma_prior;
public:
 PriorManager(const std::vector<std::string> &v_bands, const std::vector<Aerosol_Class> &v_aerosol_classes,
     const int nb_days_acc_period, const std::map<std::string, std::vector<double>> &map_prior_rpv,
     const double sigma_prior_rpv, const std::vector<double> &v_prior_AOT, const double sigma_AOT,
     const double coeff_fine_coarse_mode, const double dust_prior, const double dust_sigma_prior,
     const std::string aot_prior_directory)
   : nb_bands(v_bands.size()), v_bands(v_bands), nb_aerosol_classes(v_aerosol_classes.size()),
     v_aerosol_classes(v_aerosol_classes), sigma_AOT(sigma_AOT), dust_prior(dust_prior),
     dust_sigma_prior(dust_sigma_prior)
 {
  set_nb_days_acc_period(nb_days_acc_period);
  default_RPV_values = Eigen::MatrixXd_RowMaj::Constant(nb_bands, NB_SURFACE_PARAMS, MISSING_VALUE);
  for(int i_band = 0; i_band < nb_bands; i_band++)
  {
   default_RPV_values(i_band, RHO_0) = map_prior_rpv.at("rho_0")[i_band];
   default_RPV_values(i_band, K) = map_prior_rpv.at("k")[i_band];
   default_RPV_values(i_band, THETA) = map_prior_rpv.at("theta")[i_band];
            default_RPV_values(i_band, RHO_C) = map_prior_rpv.at("rho_c")[i_band];
  }
  default_RPV_sigma_values.setConstant(sigma_prior_rpv);
        default_RPV_sigma_values(RHO_C) = 1;
  aerosol_class_weights = Eigen::VectorXd::Ones(nb_aerosol_classes);
  int nb_aerosol_classes_fine = 0;
  int nb_aerosol_classes_coarse = 0;
  for(int i_aer = 0; i_aer < nb_aerosol_classes; i_aer++)
  {
   if(v_aerosol_classes[i_aer].type == AEROSOL_TYPE::COARSE_MODE) {
    aerosol_class_weights[i_aer] = coeff_fine_coarse_mode;
    nb_aerosol_classes_coarse++;
    }
   else nb_aerosol_classes_fine++;
  }
  double sum_weights = aerosol_class_weights.sum();
  prior_AOT = Eigen::MatrixXd_RowMaj::Constant(nb_bands, nb_aerosol_classes, 0.5);
  for(int i_band = 0; i_band < nb_bands; i_band++)
  {
   for(int i_aer = 0; i_aer < nb_aerosol_classes; i_aer++) {
    if(v_aerosol_classes[i_aer].type == AEROSOL_TYPE::COARSE_MODE) prior_AOT(i_band, i_aer) = coeff_fine_coarse_mode;
    else prior_AOT(i_band, i_aer) = (v_prior_AOT[i_band]-coeff_fine_coarse_mode)/nb_aerosol_classes_fine;
   }
  }
  prior_AOT_coarse = Eigen::MatrixXd_RowMaj::Constant(180, 360, 0.0);
  prior_AOT_fine_total = Eigen::MatrixXd_RowMaj::Constant(180, 360, 0.0);
  prior_AOT_fine_abs = Eigen::MatrixXd_RowMaj::Constant(180, 360, 0.0);
  std::unique_ptr<netCDF::NcFile> aotNCFile(new netCDF::NcFile(aot_prior_directory + "gt_c_00550nm.nc", netCDF::NcFile::read));
     aotNCFile->getVar("lat").getVar(latitudes);
     aotNCFile->getVar("lon").getVar(longitudes);
     aotNCFile->getVar("aod_ann").getVar(prior_AOT_coarse.data());
  aotNCFile.reset(new netCDF::NcFile(aot_prior_directory + "gt_f_00550nm.nc", netCDF::NcFile::read));
     aotNCFile->getVar("aod_ann").getVar(prior_AOT_fine_total.data());
  aotNCFile.reset(new netCDF::NcFile(aot_prior_directory + "gt_fa_0550nm.nc", netCDF::NcFile::read));
     aotNCFile->getVar("aer_data_ann").getVar(prior_AOT_fine_abs.data());
  #ifdef DEBUG_PRIOR
   cout << "i_band, RPV param, default value:" << endl;
   for(int i_band = 0; i_band < nb_bands; i_band++)
   {
    cout << i_band << " " << "rho_0" << " " << default_RPV_values(i_band, RHO_0) << endl;
    cout << i_band << " " << "k" << " " << default_RPV_values(i_band, K) << endl;
    cout << i_band << " " << "theta" << " " << default_RPV_values(i_band, THETA) << endl;
    cout << i_band << " " << "rho_c" << " " << default_RPV_values(i_band, RHO_C) << endl;
    cout << i_band << " " << "sigma" << " " << default_RPV_sigma_values(i_band) << endl;
    cout << endl;
   }
   cout << endl;
   cout << "i_band, i_aerosol, prior AOT:" << endl;
   for(int i_band = 0; i_band < nb_bands; i_band++)
   {
    for(int i_aer = 0; i_aer < nb_aerosol_classes; i_aer++)
    {
     cout << i_band << " " << i_aer << " " << prior_AOT(i_band, i_aer) << endl;
    }
    cout << endl;
   }
   cout << endl;
   cout << "sigma prior AOT " << sigma_AOT << endl;
   cout << endl;
  #endif
 }
 virtual ~PriorManager() { }
 void set_nb_days_acc_period(const int nb_days_acc_period)
 {
  this->nb_days_acc_period = nb_days_acc_period;
  default_increase_sigma_factor = pow(default_increase_sigma_factor_per_day, nb_days_acc_period);
 }
 void get_prior_range(const double previous_prior, const int nb_days_since_last_update, double& min_prior, double& max_prior)
 {
  min_prior = previous_prior - (max_daily_change * nb_days_since_last_update);
  max_prior = previous_prior + (max_daily_change * nb_days_since_last_update);
 }
 void get_sigma_prior_range(const int param_index, const double new_prior, double& min_sigma_prior, double& max_sigma_prior)
 {
  min_sigma_prior = factor_sigma * fabs(new_prior);
  if(min_sigma_prior < min_threshold_sigma) min_sigma_prior = min_threshold_sigma;
  max_sigma_prior = default_RPV_sigma_values(param_index);
 }
    PriorData create_default_prior(const int nb_lines, const int nb_columns, const double start_date, const std::vector<double>& latitude, const std::vector<double>& longitude, const std::string cisar_lut_directory, const std::string cisar_lut_ID)
 {
  PriorData prior_data = PriorData(nb_lines, nb_columns, v_bands, v_aerosol_classes);
  prior_data.last_update_time.setConstant(start_date);
  float ext_coeff_f_abs[prior_data.nb_bands];
  float ext_coeff_f_non_abs[prior_data.nb_bands];
  float ext_coeff_c[prior_data.nb_bands];
     float ext_coeff55_c;
     float ext_coeff55_f_abs;
     float ext_coeff55_f_non_abs;
     double da_ext055[MAX_NB_AEROSOLS_CLASSES];
     double da_ext[MAX_NB_BANDS][MAX_NB_AEROSOLS_CLASSES];
     double da_sca[MAX_NB_BANDS][MAX_NB_AEROSOLS_CLASSES];
     char ca_aero_model_name[MAX_NB_AEROSOLS_CLASSES][64];
     char ca_aero_model_type[MAX_NB_AEROSOLS_CLASSES];
     c_get_aerosol_coef(&da_ext055[0], &da_ext[0][0], &da_sca[0][0], &ca_aero_model_name[0], &ca_aero_model_type[0]);
  int position = -1;
  int position_f_abs = -1;
  int position_f_nabs = -1;
  double albedo_f_abs = 1.0;
  double albedo_f_nabs = 1.0;
  double w;
  for (Aerosol_Class &aer_class : v_aerosol_classes)
  {
   c_get_aerosol_index(aer_class.name.c_str(), &position);
   if(aer_class.type == AEROSOL_TYPE::COARSE_MODE) {
    ext_coeff55_c = da_ext055[position];
    for(int i_band=0;i_band<prior_data.nb_bands; i_band++){
     ext_coeff_c[i_band] = da_ext[i_band][position];
    }
   }
   else {
    w = da_sca[0][position]/da_ext[0][position];
    if(position_f_nabs==-1 & position_f_abs==-1){
     albedo_f_nabs = w;
     position_f_nabs = position;
     albedo_f_abs = w;
     position_f_abs = position;
    }
    else if(w<albedo_f_abs) {
     albedo_f_abs = w;
        position_f_abs = position;
    }
    else{
     albedo_f_nabs = w;
     position_f_nabs = position;
    }
   }
  }
  ext_coeff55_f_abs = da_ext055[position_f_abs];
  ext_coeff55_f_non_abs = da_ext055[position_f_nabs];
  for(int i_band=0;i_band<prior_data.nb_bands; i_band++){
   ext_coeff_f_non_abs[i_band] = da_ext[i_band][position_f_nabs];
   ext_coeff_f_abs[i_band] = da_ext[i_band][position_f_abs];
  }
  for(int i_pixel = 0; i_pixel < prior_data.nb_pixels; i_pixel++)
  {
            if (latitude[i_pixel] != -999 && longitude[i_pixel] != -999){
                int latIndex = 0;
                int lonIndex = 0;
                while(latitudes[latIndex] > latitude[i_pixel])
                {
                    latIndex ++;
                }
                while(longitudes[lonIndex] < longitude[i_pixel])
                {
                    lonIndex ++;
                }
                for(int i_band = 0 ; i_band < prior_data.nb_bands; i_band++)
                {
                    prior_data.rho_0(i_band, i_pixel) = default_RPV_values(i_band, RHO_0);
                    prior_data.k(i_band, i_pixel) = default_RPV_values(i_band, K);
                    prior_data.theta(i_band, i_pixel) = default_RPV_values(i_band, THETA);
                    prior_data.rho_c(i_band, i_pixel) = default_RPV_values(i_band, RHO_C);
                    prior_data.sigma_rho_0(i_band, i_pixel) = default_RPV_sigma_values(RHO_0);
                    prior_data.sigma_k(i_band, i_pixel) = default_RPV_sigma_values(K);
                    prior_data.sigma_theta(i_band, i_pixel) = default_RPV_sigma_values(THETA);
                    prior_data.sigma_rho_c(i_band, i_pixel) = default_RPV_sigma_values(RHO_C);
                    prior_data.aot_fine_abs(i_band,i_pixel) = prior_AOT_fine_abs(latIndex, lonIndex)/ext_coeff55_f_abs*ext_coeff_f_abs[i_band];
                    prior_data.aot_fine_non_abs(i_band,i_pixel) = (prior_AOT_fine_total(latIndex, lonIndex)-prior_AOT_fine_abs(latIndex, lonIndex))/ext_coeff55_f_non_abs*ext_coeff_f_non_abs[i_band];
                    prior_data.aot_coarse(i_band,i_pixel) = prior_AOT_coarse(latIndex, lonIndex);
                    prior_data.sigma_aot(i_band, i_pixel) = sigma_AOT;
               }
            }
        }
        return prior_data;
 }
 void set_prior_AOT(Pixel_Input& px_input, const std::vector<bool>& v_dust_pixel)
 {
  for(int i_band = 0 ; i_band < nb_bands; i_band++)
  {
   for(int i_aer = 0; i_aer < nb_aerosol_classes; i_aer++)
   {
    for(int i_obs = 0; i_obs < px_input.nb_obs; i_obs++)
    {
     if(v_dust_pixel[i_obs] && v_aerosol_classes[i_aer].type == AEROSOL_TYPE::COARSE_MODE)
     {
      px_input.prior_tau[i_band][i_aer][i_obs] = dust_prior;
      px_input.sigma_prior_tau[i_band][i_aer][i_obs] = dust_sigma_prior;
     }
     else
     {
      px_input.prior_tau[i_band][i_aer][i_obs] = prior_AOT(i_band, i_aer);
      px_input.sigma_prior_tau[i_band][i_aer][i_obs] = sigma_AOT;
     }
    }
   }
  }
 }
 void update_prior(std::vector<Pixel_Input>& v_pixel_input, std::vector<Pixel_Output>& v_pixel_output, const double end_acc_period_date, PriorData& prior_data)
 {
  if( prior_data.nb_pixels != v_pixel_input.size() || prior_data.nb_pixels != v_pixel_output.size() )
   throw std::runtime_error("PriorManager::update_prior -- nb_pixels, v_pixel_input.size(), v_pixel_output.size() and v_last_update_date.size() must be equal. Actual values are : " +
     std::to_string(prior_data.nb_pixels) + " " + std::to_string(v_pixel_input.size()) + " " + std::to_string(v_pixel_output.size()) );
  for(int i_pixel = 0; i_pixel < prior_data.nb_pixels; i_pixel++)
  {
   update_one_pixel(i_pixel, v_pixel_input[i_pixel], v_pixel_output[i_pixel], end_acc_period_date, prior_data);
  }
  prior_data.nb_days_since_last_reset += nb_days_acc_period;
        if(prior_data.nb_days_since_last_reset > prior_data.nb_days_reset)
        {
            prior_data.nb_days_since_last_reset = 0;
            prior_data.step_mean++;
            int i_reset = ( pow(-1, prior_data.step_mean+1) + 1 ) / 2;
            prior_data.nb_data.row(i_reset).setZero();
        }
 }
 void set_prior_to_fortran_input_structure(const int i_pixel, const float latitude, const float longitude, Pixel_Input& px_input, double& last_update_time,
            const PriorData& prior_data, const std::vector<bool>& v_dust_pixel, const std::string cisar_lut_directory, const std::string cisar_lut_ID)
 {
  if(i_pixel >= prior_data.nb_pixels)
   throw std::runtime_error("PriorManager::set_prior_to_fortran_input_structure -- i_pixel must be smaller than nb_pixels. Actual values are: " +
     std::to_string(i_pixel) + ", " + std::to_string(prior_data.nb_pixels));
  last_update_time = prior_data.last_update_time(i_pixel);
  double da_ext055[MAX_NB_AEROSOLS_CLASSES];
  double da_ext[MAX_NB_BANDS][MAX_NB_AEROSOLS_CLASSES];
  double da_sca[MAX_NB_BANDS][MAX_NB_AEROSOLS_CLASSES];
  char ca_aero_model_name[MAX_NB_AEROSOLS_CLASSES][64];
  char ca_aero_model_type[MAX_NB_AEROSOLS_CLASSES];
  c_get_aerosol_coef(&da_ext055[0], &da_ext[0][0], &da_sca[0][0], &ca_aero_model_name[0], &ca_aero_model_type[0]);
  int position = -1;
  std::string fine_absorbing_class = "";
  double albedo_f_abs = 1.0;
  double w;
  int nb_fine_mode_classes = 0;
  for (Aerosol_Class &aer_class : v_aerosol_classes)
  {
   c_get_aerosol_index(aer_class.name.c_str(), &position);
   if(aer_class.type == AEROSOL_TYPE::FINE_MODE) {
    nb_fine_mode_classes++;
    w = da_sca[0][position]/da_ext[0][position];
    if(fine_absorbing_class==""){
     albedo_f_abs = w;
     fine_absorbing_class = aer_class.name;
    }
    else if(w<albedo_f_abs) {
     albedo_f_abs = w;
     fine_absorbing_class = aer_class.name;
    }
   }
  }
  for(int i_band = 0 ; i_band < prior_data.nb_bands; i_band++)
  {
   px_input.prior_srf[i_band][RPV_Params_Order::RHO_0] = prior_data.rho_0(i_band, i_pixel);
   px_input.prior_srf[i_band][RPV_Params_Order::K] = prior_data.k(i_band, i_pixel);
   px_input.prior_srf[i_band][RPV_Params_Order::THETA] = prior_data.theta(i_band, i_pixel);
   px_input.prior_srf[i_band][RPV_Params_Order::RHO_C] = prior_data.rho_c(i_band, i_pixel);
   px_input.sigma_prior_srf[i_band][RPV_Params_Order::RHO_0] = prior_data.sigma_rho_0(i_band, i_pixel);
   px_input.sigma_prior_srf[i_band][RPV_Params_Order::K] = prior_data.sigma_k(i_band, i_pixel);
   px_input.sigma_prior_srf[i_band][RPV_Params_Order::THETA] = prior_data.sigma_theta(i_band, i_pixel);
   px_input.sigma_prior_srf[i_band][RPV_Params_Order::RHO_C] = prior_data.sigma_rho_c(i_band, i_pixel);
   for(int i_aer = 0; i_aer < px_input.nb_aerosols; i_aer++)
   {
    for(int i_obs = 0; i_obs < px_input.nb_obs; i_obs++)
    {
     if(v_dust_pixel[i_obs] && v_aerosol_classes[i_aer].type == AEROSOL_TYPE::COARSE_MODE)
     {
      px_input.prior_tau[i_band][i_aer][i_obs] = dust_prior;
      px_input.sigma_prior_tau[i_band][i_aer][i_obs] = dust_sigma_prior;
     }
     else
     {
      if(v_aerosol_classes.size()==1){
       px_input.prior_tau[i_band][i_aer][i_obs] = prior_data.aot_coarse(i_band, i_pixel) + prior_data.aot_fine_abs(i_band, i_pixel) + prior_data.aot_fine_non_abs(i_band, i_pixel);
      }
      else if(v_aerosol_classes[i_aer].type == AEROSOL_TYPE::COARSE_MODE)
      {
       if(prior_data.aot_coarse(i_band, i_pixel) < 0.2)
       {
        px_input.prior_tau[i_band][i_aer][i_obs] = prior_data.aot_coarse(i_band, i_pixel);
       }
       else
       {
        px_input.prior_tau[i_band][i_aer][i_obs] = 0.2;
       }
      }
      else if(nb_fine_mode_classes==1){
       px_input.prior_tau[i_band][i_aer][i_obs] = prior_data.aot_fine_abs(i_band, i_pixel) + prior_data.aot_fine_non_abs(i_band, i_pixel);
      }
      else if(v_aerosol_classes[i_aer].name == fine_absorbing_class)
      {
       px_input.prior_tau[i_band][i_aer][i_obs] = prior_data.aot_fine_abs(i_band, i_pixel);
      }
      else
      {
       px_input.prior_tau[i_band][i_aer][i_obs] = prior_data.aot_fine_non_abs(i_band, i_pixel);
      }
      px_input.sigma_prior_tau[i_band][i_aer][i_obs] = prior_data.sigma_aot(i_band,i_pixel);
     }
    }
   }
  }
 }
private:
 void update_one_pixel(const int i_pixel, const Pixel_Input& px_input, const Pixel_Output& px_output, const double end_acc_period_date, PriorData& prior_data);
};
#endif
