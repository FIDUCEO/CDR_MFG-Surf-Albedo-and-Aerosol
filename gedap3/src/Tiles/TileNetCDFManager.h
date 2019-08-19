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
#ifndef SRC_TILES_TILENETCDFMANAGER_H_
#define SRC_TILES_TILENETCDFMANAGER_H_
#include "FortranStructures.h"
#include "../Common/Global_header.h"
#include "../Common/Logger.h"
#include "../Common/Data.h"
class TileNetCDFManager
{
    std::shared_ptr<AbstractLogger> logger;
public:
 TileNetCDFManager(std::shared_ptr<AbstractLogger> logger) {
  this->logger = logger;
 }
 virtual ~TileNetCDFManager() = default;
    void load_static_data(const std::string file, StaticData& static_data);
 void load_input_data(const std::string file, InputData &input_data);
 void load_prior_data(const std::string file, PriorData &prior_data,
                         const std::vector<Aerosol_Class> &v_aerosol_classes);
 void save_static_data(const std::string file, StaticData &static_data);
 void save_input_data(const std::string file, InputData &input_data);
 void save_prior_data(const std::string file, PriorData &prior_data);
    void save_solution_data(const std::string file, std::vector<Pixel_Output> &v_pixel_output,
                            std::vector<Pixel_Input> &v_pixel_input, std::vector<Error> &v_error,
                            std::vector<std::vector<bool>> &v_dust, std::vector<std::vector<int>> &v_selected_angles,
                            const int nb_lines, const int nb_columns, const std::vector<std::string> &v_bands,
                            const std::vector<Aerosol_Class> &v_aerosol_classes,
                            const bool b_temporal_smoothness,
                            const bool b_spectral_constraint,
                            const std::vector<std::string> &satellite_images_filenames,
                            const std::vector<std::string> &list_auxiliary_data);
};
#endif
