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
#include "../Tile_Processor_Module/PriorManager.h"
void PriorManager::update_one_pixel(const int i_pixel, const Pixel_Input& px_input, const Pixel_Output& px_output, const double end_acc_period_date, PriorData& prior_data)
{
 int nb_days_since_last_update = (int) (end_acc_period_date - prior_data.last_update_time[i_pixel] + 0.001);
 double new_prior = 0.0;
 double previous_prior = 0.0;
 double min_prior = 0.0;
 double max_prior = 0.0;
 double new_sigma_prior = 0.0;
 double min_sigma_prior = 0.0;
 double max_sigma_prior = 0.0;
    double cosg = 0.0;
    double tolerance = cos(15.0 * COEFF_DEG2RAD);
    bool hot_spot = false;
 int i_mean = ( pow(-1, prior_data.step_mean) + 1 ) / 2;
    double average_quality = px_output.quality[0];
    if( px_input.srf_type == LAND )
    {
       for(int i_angle=0;i_angle<px_input.nb_angles;i_angle++){
            cosg = cos(px_input.SZA[i_angle]) * cos(px_input.VZA[i_angle]) + sin(px_input.SZA[i_angle]) * sin(px_input.VZA[i_angle]) * cos(px_input.RAA[i_angle]);
            if(cosg>tolerance){
                hot_spot = true;
                break;
            }
        }
        for(int i_rpv_param = 0; i_rpv_param < NB_SURFACE_PARAMS; i_rpv_param++)
        {
            for(int i_band = 0; i_band < px_input.nb_bands; i_band++)
            {
                if(px_output.status == 0) {
                    new_prior = px_output.srf_param[i_band][i_rpv_param];
                    previous_prior = px_input.prior_srf[i_band][i_rpv_param];
                    get_prior_range(previous_prior, nb_days_since_last_update, min_prior, max_prior);
                    if(new_prior < min_prior) new_prior = min_prior;
                    if(new_prior > max_prior) new_prior = max_prior;
                    if (average_quality>0.0){
                        for(int i_vector = 0; i_vector < prior_data.nb_mean_vectors; i_vector++)
                        {
                            if(prior_data.nb_data(i_vector, i_pixel) == 0)
                            {
                                prior_data.get_sum(i_vector, i_rpv_param, i_band, i_pixel) = 0.0;
                                prior_data.quality(i_vector, i_pixel) = 0.0;
                                prior_data.get_min_rpv(i_vector, i_rpv_param, i_band, i_pixel) = px_output.srf_param[i_band][i_rpv_param];
                                prior_data.get_max_rpv(i_vector, i_rpv_param, i_band, i_pixel) = px_output.srf_param[i_band][i_rpv_param];
                            }
                            prior_data.get_sum(i_vector, i_rpv_param, i_band, i_pixel) += (px_output.srf_param[i_band][i_rpv_param]*average_quality);
                            if(px_output.srf_param[i_band][i_rpv_param] > prior_data.get_max_rpv(i_vector, i_rpv_param, i_band, i_pixel)){
                                prior_data.get_max_rpv(i_vector, i_rpv_param, i_band, i_pixel) = px_output.srf_param[i_band][i_rpv_param];
                            }
                            if(px_output.srf_param[i_band][i_rpv_param] < prior_data.get_min_rpv(i_vector, i_rpv_param, i_band, i_pixel)) {
                                prior_data.get_min_rpv(i_vector, i_rpv_param, i_band, i_pixel) = px_output.srf_param[i_band][i_rpv_param];
                            }
                            prior_data.quality(i_vector, i_pixel) +=average_quality;
                       }
                    }
                    if( prior_data.nb_data(i_mean, i_pixel)+1.0 >= prior_data.min_data &&(i_rpv_param!=RPV_Params_Order::RHO_0 && i_rpv_param!=RPV_Params_Order::RHO_C))
                    {
                        new_prior = prior_data.get_sum(i_mean, i_rpv_param, i_band, i_pixel) / prior_data.quality(i_mean, i_pixel);
                    }
                    check_RPV_physical_range(i_rpv_param, new_prior);
                    switch (i_rpv_param) {
                    case RHO_0:
                        prior_data.rho_0(i_band, i_pixel) = new_prior;
                        break;
                    case K:
                        prior_data.k(i_band, i_pixel) = new_prior;
                        break;
                    case THETA:
                        prior_data.theta(i_band, i_pixel) = new_prior;
                        break;
                    case RHO_C:
                        prior_data.rho_c(i_band, i_pixel) = new_prior;
                        break;
                    default:
                        break;
                    }
                }
                if(px_output.status == 0)
                {
                    new_sigma_prior = px_output.sigma_srf_param[i_band][i_rpv_param];
                }
                else
                {
                    new_sigma_prior = px_input.sigma_prior_srf[i_band][i_rpv_param] * default_increase_sigma_factor;
                }
                get_sigma_prior_range(i_rpv_param, new_prior, min_sigma_prior, max_sigma_prior);
                if(new_sigma_prior < min_sigma_prior) new_sigma_prior = min_sigma_prior;
                if(new_sigma_prior > max_sigma_prior) new_sigma_prior = max_sigma_prior;
                if( prior_data.nb_data(i_mean, i_pixel) >= prior_data.min_data &&(i_rpv_param!=RPV_Params_Order::RHO_0))
                {
                    new_sigma_prior = (prior_data.get_max_rpv(i_mean, i_rpv_param, i_band, i_pixel) -
                                       prior_data.get_min_rpv(i_mean, i_rpv_param, i_band, i_pixel)) / 2.0;
                    if (new_sigma_prior>max_threshold_sigma) new_sigma_prior = max_threshold_sigma;
                    if (new_sigma_prior<min_threshold_sigma) new_sigma_prior = min_threshold_sigma;
                }
                if (hot_spot && i_rpv_param!=RPV_Params_Order::RHO_C) new_sigma_prior = (min_threshold_sigma*2)/cosg;
                switch (i_rpv_param) {
                case RHO_0:
                    prior_data.sigma_rho_0(i_band, i_pixel) = new_sigma_prior;
                    break;
                case K:
                    prior_data.sigma_k(i_band, i_pixel) = new_sigma_prior;
                    break;
                case THETA:
                    prior_data.sigma_theta(i_band, i_pixel) = new_sigma_prior;
                    break;
                case RHO_C:
                    prior_data.sigma_rho_c(i_band, i_pixel) = new_sigma_prior;
                    break;
                default:
                    break;
                }
                if(px_output.status == 0)
                    prior_data.last_update_time[i_pixel] = end_acc_period_date;
            }
        }
        if(average_quality>0.0){
            for(int i_vector = 0; i_vector < prior_data.nb_mean_vectors; i_vector++)
                prior_data.nb_data(i_vector, i_pixel)++;
        }
    }
 else
 {
  for(int i_band = 0; i_band < px_input.nb_bands; i_band++)
  {
   prior_data.rho_0(i_band, i_pixel) = MISSING_VALUE;
   prior_data.sigma_rho_0(i_band, i_pixel) = 0.0;
   prior_data.k(i_band, i_pixel) = MISSING_VALUE;
   prior_data.sigma_k(i_band, i_pixel) = 0.0;
   prior_data.theta(i_band, i_pixel) = MISSING_VALUE;
   prior_data.sigma_theta(i_band, i_pixel) = 0.0;
   prior_data.rho_c(i_band, i_pixel) = MISSING_VALUE;
   prior_data.sigma_rho_c(i_band, i_pixel) = 0.0;
  }
 }
}
