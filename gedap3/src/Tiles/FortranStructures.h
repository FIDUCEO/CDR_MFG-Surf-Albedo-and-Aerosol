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
#ifndef SOURCE_TILE_PROCESSOR_MODULE_FORTRANSTRUCTURES_H_
#define SOURCE_TILE_PROCESSOR_MODULE_FORTRANSTRUCTURES_H_
#include "../Common/Global_header.h"
using std::cout;
using std::endl;
 typedef double real_t;
enum Surface_Type { LAND = 1, SEA };
struct Pixel_Input
{
 int max_nb_iterations = 10;
 int nb_bands = 0;
 int nb_aerosols = 0;
 int nb_angles = 0;
 int nb_obs = 0;
 int srf_type = -1;
 real_t SZA[MAX_NB_ANGLES];
 real_t VZA[MAX_NB_ANGLES];
 real_t RAA[MAX_NB_ANGLES];
 real_t time[MAX_NB_ANGLES];
 real_t SrfPres[MAX_NB_ATM];
 real_t WV_COL[MAX_NB_ATM];
 real_t OZ_COL[MAX_NB_ATM];
 real_t layer_height[MAX_NB_ATM];
 real_t wind_speed = 0.0;
 real_t wind_direction = 0.0;
 int aerosol_mask[MAX_NB_AEROSOLS_CLASSES];
 real_t TOA_BRF[MAX_NB_BANDS][MAX_NB_ANGLES];
 real_t sigma_TOA_BRF[MAX_NB_BANDS][MAX_NB_ANGLES];
    real_t cloud_mask[MAX_NB_ANGLES];
 real_t first_guess_srf[MAX_NB_BANDS][NB_SURFACE_PARAMS];
 real_t first_guess_tau[MAX_NB_BANDS][MAX_NB_AEROSOLS_CLASSES][MAX_NB_ANGLES];
 real_t prior_srf[MAX_NB_BANDS][NB_SURFACE_PARAMS];
 real_t sigma_prior_srf[MAX_NB_BANDS][NB_SURFACE_PARAMS];
 real_t prior_tau[MAX_NB_BANDS][MAX_NB_AEROSOLS_CLASSES][MAX_NB_ANGLES];
 real_t sigma_prior_tau[MAX_NB_BANDS][MAX_NB_AEROSOLS_CLASSES][MAX_NB_ANGLES];
 real_t sigma_temp_smoothness[MAX_NB_BANDS][MAX_NB_ANGLES];
 real_t sigma_spectral_constraint = 0.05;
 real_t sigma_surface_spectral_constraint[MAX_NB_BANDS][NB_SURFACE_PARAMS];
 Pixel_Input()
 {
  for(int i_angle = 0; i_angle < MAX_NB_ANGLES; i_angle++)
   time[i_angle] = MISSING_VALUE;
  for(int i_band = 0; i_band < MAX_NB_BANDS; i_band++)
   for(int i_angle = 0; i_angle < MAX_NB_ANGLES; i_angle++)
    sigma_temp_smoothness[i_band][i_angle] = 0.05;
 }
 Pixel_Input(const int nb_bands, const int nb_aerosols, const int nb_angles, const int nb_obs) :
  nb_bands(nb_bands), nb_aerosols(nb_aerosols), nb_angles(nb_angles), nb_obs(nb_obs)
 {
  for(int i_angle = 0; i_angle < MAX_NB_ANGLES; i_angle++)
   time[i_angle] = MISSING_VALUE;
  for(int i_band = 0; i_band < MAX_NB_BANDS; i_band++)
   for(int i_angle = 0; i_angle < nb_obs; i_angle++)
    sigma_temp_smoothness[i_band][i_angle] = 0.05;
 }
 ~Pixel_Input() { }
    void print() const
 {
 }
};
struct Pixel_Output
{
 int nb_bands = 0;
 int nb_aerosols = 0;
 int nb_angles = 0;
 int nb_obs = 0;
 int nb_iterations = -1;
 real_t quality[MAX_NB_ANGLES];
 real_t quality_information[NB_QUALITY_TEST][MAX_NB_ANGLES];
 int status = -1;
 real_t time[MAX_NB_ANGLES] ;
 real_t cost = MISSING_VALUE;
 real_t srf_param[MAX_NB_BANDS][NB_SURFACE_PARAMS];
 real_t sigma_srf_param[MAX_NB_BANDS][NB_SURFACE_PARAMS];
 real_t tau[MAX_NB_BANDS][MAX_NB_AEROSOLS_CLASSES][MAX_NB_ANGLES];
 real_t sigma_tau[MAX_NB_BANDS][MAX_NB_AEROSOLS_CLASSES][MAX_NB_ANGLES];
 real_t total_tau[MAX_NB_BANDS][MAX_NB_ANGLES];
 real_t sigma_total_tau[MAX_NB_BANDS][MAX_NB_ANGLES];
 real_t tau_55[MAX_NB_AEROSOLS_CLASSES][MAX_NB_ANGLES];
 real_t sigma_tau_55[MAX_NB_AEROSOLS_CLASSES][MAX_NB_ANGLES];
 real_t SSA[MAX_NB_BANDS][MAX_NB_ANGLES];
 real_t asym[MAX_NB_BANDS][MAX_NB_ANGLES];
 real_t sigma_SSA[MAX_NB_BANDS][MAX_NB_ANGLES];
 real_t sigma_asym[MAX_NB_BANDS][MAX_NB_ANGLES];
 real_t TOA_BRF[MAX_NB_BANDS][MAX_NB_ANGLES];
 real_t BHR[MAX_NB_BANDS];
 real_t sigma_BHR[MAX_NB_BANDS];
    Pixel_Output()
 {
  for(int i_angle = 0; i_angle < MAX_NB_ANGLES; i_angle++)
   time[i_angle] = MISSING_VALUE;
 }
 Pixel_Output(const int nb_bands, const int nb_aerosols, const int nb_angles, const int nb_obs) :
  nb_bands(nb_bands), nb_aerosols(nb_aerosols), nb_angles(nb_angles), nb_obs(nb_obs)
 {
  for(int i_angle = 0; i_angle < MAX_NB_ANGLES; i_angle++)
   time[i_angle] = MISSING_VALUE;
 }
 void print()
 {
 }
};
#endif /* SOURCE_TILE_PROCESSOR_MODULE_FORTRANSTRUCTURES_H_ */
