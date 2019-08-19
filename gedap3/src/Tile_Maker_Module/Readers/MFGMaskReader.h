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
 * MVIRIMaskReader.h
 *
 *  Created on: 27 Jul 2017
 *      Author: cedric
 */

#ifndef TILE_MAKER_MODULE_MFGMASKREADER_H_
#define TILE_MAKER_MODULE_MFGMASKREADER_H_

#include <netcdf>

#include "AbstractReader.h"
#include "Satellite_Parameters.h"

//#define DEBUG_MASK


/**
 * \ingroup Tile_Maker
 *
 * \brief   Read external masks and aggregate them in one unique mask. The data are saved in an instance of Mask structure.
 * \details This class use the STL bitset class. \sa Mask
 */
class MFGMaskReader : public AbstractReader
{
    // if IODC the file name doesn't change
    std::string m_path_file_clouds;

	std::shared_ptr<netCDF::NcFile> inputFile;
    int m_reliability_cfc_retrieval;

    Eigen::ArrayXXsc_RowMaj pcfc;
    Eigen::ArrayXXi_RowMaj cfc;
    Eigen::ArrayXXs_RowMaj cmaskscore;

public:
    MFGMaskReader(std::shared_ptr<AbstractLogger> logger);

    ~MFGMaskReader() override = default;

    /**
     * @brief setup
     * @param inputFilename
     * @param logger
     */

    void setup(const std::string inputFilename, const int cloud_free_probability,
			   std::shared_ptr<AbstractLogger> logger);

    /**
	 * Reads external mask files and store the data
	 */
	void read();

    void IODC_read(const int slot);

	/**
	 * For a given tile, aggregates the external masks into one and copy it into an instance of Mask structure
	 */
	void extract_data(Mask& mask_data, const bool useCMASKSCORE, const bool usePCFC, TileDimensions& tile_dims);

};

#endif /* TILE_MAKER_MODULE_MFGMASKREADER_H_ */
