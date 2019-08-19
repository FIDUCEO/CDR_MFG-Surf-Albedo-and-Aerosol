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
 * MVIRIMaskReader.cpp
 *
 *  Created on: 27 Jul 2017
 *      Author: Cedric
 */

#include "MFGMaskReader.h"

using namespace netCDF;
using exceptions::NcException;

MFGMaskReader::MFGMaskReader(std::shared_ptr<AbstractLogger> logger)
        : AbstractReader(logger),
          cfc(IMAGE_NB_LINES, IMAGE_NB_COLUMNS),
          pcfc(IMAGE_NB_LINES, IMAGE_NB_COLUMNS),
          cmaskscore(IMAGE_NB_LINES, IMAGE_NB_COLUMNS) {

    // Initialize buffers with zeros
    cfc = Eigen::ArrayXXi_RowMaj::Zero(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
    pcfc = Eigen::ArrayXXsc_RowMaj::Zero(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
    cmaskscore = Eigen::ArrayXXs_RowMaj::Zero(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);

}

void MFGMaskReader::setup(const std::string inputFilename, const int reliability_cfc_retrieval, std::shared_ptr<AbstractLogger> logger) {

    // Try to open input file
    isReady = true;
    if (inputFilename.compare(m_path_file_clouds) != 0)
    {
        try {
            logger->write_log(LogLevel::_INFO_, "MFGMaskReader", "setup",
                              "Opening mask file " + inputFilename);
            inputFile = std::make_shared<NcFile>(inputFilename, NcFile::read);

            m_path_file_clouds = inputFilename;
        }
        catch (NcException &c) {

            logger->write_log(LogLevel::_ERROR_, "MFGMaskReader", "setup",
                              "Error while opening cloud mask file " + inputFilename);
            isReady = false;
        }

        m_reliability_cfc_retrieval = reliability_cfc_retrieval;

        // Signal that cache was successfuly updated
        logger->write_log(LogLevel::_DEBUG_, "MFGMaskReader", "setup",
                          "Cloud mask opened.");

    }
    else {
        logger->write_log(LogLevel::_DEBUG_, "MFGMaskReader", "setup",
                          "Nothing to do, IODC daily mask loaded");
    }

}


void MFGMaskReader::read() {

    Eigen::ArrayXXs_RowMaj packed_values_short(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
    Eigen::ArrayXXi_RowMaj packed_values_signed_integer_cfc(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
    Eigen::ArrayXXi_RowMaj packed_values_signed_integer_pcfc(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);

    float add_offset = 0.;
    float scale_factor = 0.;

    try {
        // read and unpack cfc, pcfc and cmaskscore
        inputFile->getVar("CFC").getVar(packed_values_signed_integer_cfc.data());
        for (int line = 0; line < IMAGE_NB_LINES; line++) {
            for (int column = 0; column < IMAGE_NB_COLUMNS; column++) {
                cfc(line, column) = (packed_values_signed_integer_cfc(line, column));
            }
        }

        inputFile->getVar("PCFC").getVar(packed_values_signed_integer_pcfc.data());
        for (int line = 0; line < IMAGE_NB_LINES; line++) {
            for (int column = 0; column < IMAGE_NB_COLUMNS; column++) {
                pcfc(line, column) = (packed_values_signed_integer_pcfc(line, column));
            }
        }

        // get values of scale factor and offset
        inputFile->getVar("CMASKSCORE").getVar(packed_values_short.data());
        inputFile->getVar("CMASKSCORE").getAtt("add_offset").getValues(&add_offset);
        inputFile->getVar("CMASKSCORE").getAtt("scale_factor").getValues(&scale_factor);

        for (int line = 0; line < IMAGE_NB_LINES; line++) {
            for (int column = 0; column < IMAGE_NB_COLUMNS; column++) {
                cmaskscore(line, column) = (packed_values_short(line, column) * scale_factor) + add_offset;
            }
        }

    }
    catch (NcException &c) {
        logger->write_log(LogLevel::_ERROR_, "MFGMaskReader", "read",
                          "error when reading cloud mask file");
    }

}



//***************************************************************************
// prepare a specific reader for the IODC

void MFGMaskReader::IODC_read(const int slot)
{

        // Allocate buffers
        std::vector<size_t> start = {static_cast<size_t>(slot),
                                     0,
                                     0};
        std::vector<size_t> count = {1,
                                     static_cast<size_t>(IMAGE_NB_LINES),
                                     static_cast<size_t>(IMAGE_NB_COLUMNS)};


        Eigen::ArrayXXs_RowMaj packed_values_short(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
        Eigen::ArrayXXi_RowMaj packed_values_signed_integer_cfc(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);
        Eigen::ArrayXXi_RowMaj packed_values_signed_integer_pcfc(IMAGE_NB_LINES, IMAGE_NB_COLUMNS);

        float add_offset = 0.;
        float scale_factor = 0.;

        try {
            // read and unpack cfc, pcfc and cmaskscore

            inputFile->getVar("CFC").getVar(start, count, packed_values_signed_integer_cfc.data());

            for (int line = 0; line < IMAGE_NB_LINES; line++){
                for (int column = 0; column < IMAGE_NB_COLUMNS; column++){
                    cfc(line, column) = (packed_values_signed_integer_cfc(line, column));
                }
            }

            inputFile->getVar("PCFC").getVar(start, count, packed_values_signed_integer_pcfc.data());
            for (int line = 0; line < IMAGE_NB_LINES; line++) {
                for (int column = 0; column < IMAGE_NB_COLUMNS; column++) {
                    pcfc(line, column) = (packed_values_signed_integer_pcfc(line, column));
                }
            }

            // get values of scale factor and offset
            inputFile->getVar("CMASKSCORE").getVar(start, count, packed_values_short.data());
            inputFile->getVar("CMASKSCORE").getAtt("scale_factor").getValues(&scale_factor);
            inputFile->getVar("CMASKSCORE").getAtt("add_offset").getValues(&add_offset);

            for (int line = 0; line < IMAGE_NB_LINES; line++) {
                for (int column = 0; column < IMAGE_NB_COLUMNS; column++) {
                    cmaskscore(line, column) = (packed_values_short(line, column) * scale_factor) + add_offset;
                }
            }

        }
        catch (NcException &c) {
            logger->write_log(LogLevel::_ERROR_, "MFGMaskReader", "read",
                              "error when reading cloud mask file");
        }

}


//***************************************************************************


void MFGMaskReader::extract_data(Mask &mask_data, const bool useCMASKSCORE, const bool usePCFC,
                                   TileDimensions &tile_dims) {

    /* ****************************	*
     * 		BUILD TILE CLOUD MASK			*
     * ****************************	*/
    bool isCloud;
    int i_image_line, i_image_column;

    for (int i_line = 0; i_line < mask_data.nb_lines; i_line++) {
        for (int i_column = 0; i_column < mask_data.nb_columns; i_column++) {
            i_image_line = tile_dims.first_line + i_line;
            i_image_column = tile_dims.first_column + i_column;

            // For any cfc value different than 0, the pixel is cloudy.
            // The external disk mask is set to -127, and the cfc field is a percentage ranging from 0 to 100.
            // Only pixels with cfc=0 are considered truly cloud free.
            isCloud = (cfc(i_image_line, i_image_column) != 0);

            // if the field is not a cloud and the user has defined to also use the cmaskscore to check for cloud presence,
            // the pixel will be a cloud if this value is larger than zero
            if (!isCloud and useCMASKSCORE) {
                isCloud = (cmaskscore(i_image_line, i_image_column) > 0.0);
            }

            // if the field is not a cloud and the user wants to use the field "pcfc" to check for quality of the retrival
            // the higher the value of the "pfcf" term, the more reliable the "cfc". Therefore we select a "pcfc" of 20%
            // and state that if the pixel results cloud free but the probability of actually being cloud free is low
            // we set this pixel as being cloudy.
            if (!isCloud and usePCFC) {
                isCloud = (pcfc(i_image_line, i_image_column) < m_reliability_cfc_retrieval);
            }

            std::bitset<Mask::NB_MASK_VALUES> maskValue;
            if (isCloud) maskValue.set(Mask::CLOUDS);

            mask_data.mask(i_line, i_column) = Mask::bitset_to_uchar(maskValue);
        }
    }
}





