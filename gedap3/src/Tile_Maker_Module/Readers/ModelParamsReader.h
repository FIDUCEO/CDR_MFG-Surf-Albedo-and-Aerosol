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
 * ModelParamsReader.h
 *
 *  Created on: 13 Oct 2015
 *      Author: alix
 */

#ifndef TILE_MAKER_MODULE_MODELPARAMSREADER_H_
#define TILE_MAKER_MODULE_MODELPARAMSREADER_H_

#include <netcdf>

#include "AbstractReader.h"
#include "Satellite_Parameters.h"

//#define DEBUG_MODEL_PARAMS
//#define DEBUG_MODEL_PARAMS_INTERPOLATION

/**
 * \ingroup Tile_Maker
 *
 * \brief   Structure used to interpolate ECMWF data over time (linear interpolation).
 * \details Example (SEVIRI):\n
 *          To get the value of the surface pressure for the datetime 2008-01-01 14:30, we:\n
 *          1) read the surface pressure from the 2 ECMWF files at datetime 2008-01-01 12:00 and 2008-01-01 18:00\n
 *          2) compute the gradient:
 *             \f[ \nabla_{srfpr} = \frac{srfPr(18:00) - srfPr(12:00)}{24} \f]
 *          3) In class ModelParamsReader, compute:
 *             \f[ srfPr(14:30) = srfPr(12:00) + nbsteps(12:00,14:30) * \nabla_{srfpr} \f]
 *             where \f$ nbsteps(12:00,14:30) = 2h:30 / 15 = 150 / 15 = 10 \f$.
 *
 */
struct ModelParamsGradient
{
	//! number of SEVIRI images (slots) between 2 ECMWF files (One ECMWF file every 6 hour & One SEVIRI image every 15 minutes --> (6*60)/15 = 24)
    //! number of MVIRI images (slots) between 2 ECMWF files (One ECMWF file every 6 hour & One MVIRI image every 30 minutes --> (6*60)/30 = 12)

    const int nb_steps = 12;

	Eigen::MatrixXf_RowMaj currentData;
	Eigen::MatrixXf_RowMaj nextData;
	Eigen::MatrixXf_RowMaj gradient;
    Eigen::MatrixXf_RowMaj stepData;

	ModelParamsGradient()
	{
        currentData = Eigen::MatrixXf_RowMaj(ERA_INTERIM_IMAGE_NB_LINES, ERA_INTERIM_IMAGE_NB_COLUMNS);
        nextData    = Eigen::MatrixXf_RowMaj(ERA_INTERIM_IMAGE_NB_LINES, ERA_INTERIM_IMAGE_NB_COLUMNS);
        gradient    = Eigen::MatrixXf_RowMaj(ERA_INTERIM_IMAGE_NB_LINES, ERA_INTERIM_IMAGE_NB_COLUMNS);
        stepData    = Eigen::MatrixXf_RowMaj(ERA_INTERIM_IMAGE_NB_LINES, ERA_INTERIM_IMAGE_NB_COLUMNS);
	}

	~ModelParamsGradient() = default;

	void reset_gradient() { gradient.setZero(); }

    void set_next_data_to_current()
    {
        currentData = nextData;
        nextData.setZero();
    }

// ERA-INTERIM images have a time frequence of 6 hours. For MVIRI this means that there are 12 MVIRI images between two ERA-INTERIM ones
// the gradient is the delta value of data inbetween these frames. Therefore to obtain the value at the time slot for the MVIRI image
// you neet to sum to the ERAinterim data the corresponding delta, obtained by multipling the gradient with the slot of the image you are considering
// the steps  reset_gradient,  set_next_data_to_current are compute_and_set_gradient are required only at the beginning and in

	void compute_and_set_gradient() { gradient = ( nextData - currentData ) / nb_steps; }

    void evaluate_data_at_step(const int slot) { stepData = currentData + ( slot * gradient ); }

	void print(const std::string name, int first_line, int last_line, int first_column, int last_column)
	{
		int nb_lines   = last_line - first_line + 1;
		int nb_columns = last_column - first_column + 1;

		std::cout << "Model Parameters - " << name << std::endl;

		std::cout << "Current Data: " << std::endl;
		std::cout << currentData.block(first_line, first_column, nb_lines, nb_columns) << " ";
		std::cout << std::endl;

		std::cout << "Next Data: " << std::endl;
		std::cout << nextData.block(first_line, first_column, nb_lines, nb_columns) << " ";
		std::cout << std::endl;

		std::cout << "Gradient: " << std::endl;
		std::cout << gradient.block(first_line, first_column, nb_lines, nb_columns) << " ";
		std::cout << std::endl;
	}
};


/**
 * \ingroup Tile_Maker
 *
 * \brief   Read the total column ozone (TCO3), total column water vapor (TCWV), surface pressure, wind speed and wind direction
 *          from ECMWF ERA-INTERIM files and copy the data in an instance of ModelParams structure.
 * \details \sa ModelParams
 *
 * \warning External data are assumed to be already in the right projection.\n
 *          The ECMWF data (TCO3, TCWV, surface pressure, wind speed, wind direction) are projected on the SEVIRI grid using the program SOSIE located in:\n
 *          svn://rayser01/RAYFERENCE/Algorithms/GEDAP/SOSIE_ECMWF_interpolation
 *
 * \todo    Strongly simplify this class by processing the time interpolation of ECMWF ERA-INTERIM files before to run GEDAP
 *          (e.g. for SEVIRI: create a ECMWF file for every 15 minutes)
 */
class ModelParamsReader : public AbstractReader
{
private:
    // input ERA-INTERIM netCDF file
    std::shared_ptr<netCDF::NcFile> mp_inputFile;
    Eigen::ArrayXf_RowMaj m_era_interim_latitude_cache;
    Eigen::ArrayXf_RowMaj m_era_interim_longitude_cache;

    // Static data filenames
    std::string m_geopotential_filename;
    std::string m_elevation_filename;
    std::string m_aerosol_height_filename_1;
    std::string m_aerosol_height_filename_2;

    // Cache ERA-INTERIM
    Eigen::ArrayXXs_RowMaj packed_values;
    Eigen::ArrayXXs_RowMaj packed_values_bis;

    // Static data cache
    Eigen::ArrayXXf_RowMaj m_latitude_cache;
    Eigen::ArrayXXf_RowMaj m_longitude_cache;
    Eigen::ArrayXXf_RowMaj m_geopotential_cache;
    Eigen::ArrayXXf_RowMaj m_elevation_cache;
    Eigen::ArrayXXf_RowMaj m_aerosol_height_cache_1;
    Eigen::ArrayXXf_RowMaj m_aerosol_height_cache_2;

    // Cache flag
    bool m_cache_ok;

	ModelParamsGradient TCO3;
	ModelParamsGradient TCWV;
	ModelParamsGradient SrfPress;
	ModelParamsGradient WindSpeed;
	ModelParamsGradient WindDir;
	ModelParamsGradient U10_era_image;
	ModelParamsGradient V10_era_image;

    //caching right ERA-INTERIM frame
    int slot = 0.0;

public:
    ModelParamsReader(std::shared_ptr<AbstractLogger> logger) : AbstractReader(logger), m_cache_ok(false)
    {
    }

    virtual ~ModelParamsReader() { }
//    virtual ~ModelParamsReader() {  logger->write_log(LogLevel::_DEBUG_, "ModelParamsReader", "~ModelParamsReader",
//                                                           "Calling destructor ModelParamsReader");}



	//! For testing
	ModelParamsGradient& get_wind_speed_gradiant() { return WindSpeed; }

	/**
	 * \brief compute the values of the model parameters using a linear interpolation of ECMWF ERA-INTERIM data.
	 *        Copy the data in an instance of ModelParams structure.
	 * \details \sa ModelParamsGradient
	 *
	 * \param step  number of slots since the last ECMWF file (e.g. for slot 15:00, step = (15:00 - 12:00) / 15 = 12).
	 */
    void extract_data(ModelParams &model_params_data, const std::string elevation_file,
                      const std::string geopotential_file, const std::string aerosol_height_file_1,
                      const std::string aerosol_height_file_2, const int day, const bool heightIsConstant,
                      TileDimensions &tile_dims, const int step,
                      const bool ecmwf_constant);

	/**
	 * \brief    Read the first ERA-INTERIM file of the accumulation period.
     * \warning  Must be called each time you start a new accumulation period.
	 */
    void read_first_model_params_data(const std::string inputFile, const int slot);

	/**
	 * \sa ModelParamsGradient
	 */
    void compute_gradient_for_interpolation(const std::string inputFile, const int slot);
    void create_full_era_interim_image_at_step(const int slot);

private:
	/**
	 * \note  \li Surface Pressure values are divided by 100 to get them in hPa instead of Pa.
	 *        \li Total Column Ozone concentration values are converted from \f$ kg ~ m^{-2} \f$ to DU (Dobson Unit).
	 *        \li Corrupted value in ERA-INTERIM files = -32767.
	 */
//	void read_next_data(const std::string inputFile);
    void load_era_interim_cache(const std::string inputFile, const int slot);

    /// If needed, this method loads aerosol_height, geopotential and elevation files contents into cache.
    void update_cache(const bool ecmwf_constant);
    float distance_geo(const float lat1, const float lon1, const float lat2, const float  lon2);
};

#endif /* TILE_MAKER_MODULE_MODELPARAMSREADER_H_ */
