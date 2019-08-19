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
 * CloudMaskReader.cpp
 *
 *  Created on: 13 Oct 2015
 *      Author: alix
 */

#include "MaskReader.h"

using namespace H5;

MaskReader::MaskReader(const std::string inputFile, std::shared_ptr<AbstractLogger> logger)
	: AbstractReader(logger)
{
	// Try to open input file
	isReady = true;

	this->inputFile = inputFile;

	try
	{
        HDFinputFile = std::make_shared<H5File>(inputFile, H5F_ACC_RDONLY);
	}
	catch (FileIException &error)
	{
		logger->write_log(LogLevel::_ERROR_, "MaskReader", "constructor",
						  "can't open mask file with path " + inputFile);
		isReady = false;
	}
}


MaskReader::~MaskReader()
{
	if( HDFinputFile.use_count() > 0 ) HDFinputFile->close();
}


void MaskReader::read()
{
    try{
		std::string datasetPath;
		DataSet dataset;

		/* ****************************	*
		 * 			CLOUDS				*
		 * ****************************	*/
		datasetPath = "/CMa";
		dataset = HDFinputFile->openDataSet( datasetPath );
		dataset.read(cloud_mask, PredType::NATIVE_UCHAR);

		/* ****************************	*
		 * 			DUST				*
		 * ****************************	*/
		datasetPath = "/CMa_DUST";
		dataset = HDFinputFile->openDataSet( datasetPath );
		dataset.read(dust_mask, PredType::NATIVE_UCHAR);

		/* ****************************	*
		 * 		VOLCANIC PLUME			*
		 * ****************************	*/
		datasetPath = "/CMa_VOLCANIC";
		dataset = HDFinputFile->openDataSet( datasetPath );
		dataset.read(volcanic_plume_mask, PredType::NATIVE_UCHAR);

		/* ****************************	*
		 * 			TESTS				*
		 * ****************************	*/
		datasetPath = "/CMa_TEST";
		dataset = HDFinputFile->openDataSet( datasetPath );
		dataset.read(test_mask, PredType::NATIVE_USHORT);

	}
	// catch failure caused by the DataSet operations
	catch( DataSetIException& error )
	{
	  error.printErrorStack();
	  logger->write_log(LogLevel::_ERROR_, "MaskReader", "read", "can't read mask file with path " + inputFile);
	}
	// catch failure caused by the DataSpace operations
	catch( DataSpaceIException& error )
	{
	  error.printErrorStack();
	  logger->write_log(LogLevel::_ERROR_, "MaskReader", "read", "can't read mask file with path " + inputFile);
	}
	// catch failure caused by the DataSpace operations
	catch( DataTypeIException& error )
	{
	  error.printErrorStack();
	  logger->write_log(LogLevel::_ERROR_, "MaskReader", "read", "can't read mask file with path " + inputFile);
	}
}


void MaskReader::extract_data(Mask& mask_data, TileDimensions& tile_dims)
{
	/* ****************************	*
	 * 		BUILD TILE MASK			*
	 * ****************************	*/
	bool isCloud = false;
	bool isSnow  = false;
	bool isFire  = false;

	/*
	bool isSmoke = false;
	bool testSmokeActivated = false;

	int nb_SpatialSmoothing = 0;
	for(int line = tile_dims.first_line; line <= tile_dims.last_line; line++)
	{
		for(int column = tile_dims.first_column; column <= tile_dims.last_column; column++)
		{
			if(testSmokeActivated && (cloud_mask[line][column] == 2 or cloud_mask[line][column] == 3) )
			{
				bitset<16> test_value(test_mask[line][column]);
				nb_SpatialSmoothing = bitset<16>(test_mask[line-1][column-1]).test(7) +
									  bitset<16>(test_mask[line-1][column]).test(7)   +
									  bitset<16>(test_mask[line-1][column+1]).test(7) +
									  bitset<16>(test_mask[line][column-1]).test(7) +
									  bitset<16>(test_mask[line][column+1]).test(7) +
									  bitset<16>(test_mask[line+1][column-1]).test(7) +
									  bitset<16>(test_mask[line+1][column]).test(7) +
									  bitset<16>(test_mask[line+1][column+1]).test(7);
				smoke_mask[line][column] = (test_value.test(3) or test_value.test(6) or nb_SpatialSmoothing >= 3);
			}else{
				smoke_mask[line][column] = false;
			}
		}
	}
	*/

	int i_image_line = 0;
	int i_image_column = 0;
	for(int i_line = 0; i_line < mask_data.nb_lines; i_line++)
	{
		for(int i_column = 0; i_column < mask_data.nb_columns; i_column++)
		{
			i_image_line   = tile_dims.first_line + i_line;
			i_image_column = tile_dims.first_column + i_column;

			isCloud = ((cloud_mask[i_image_line][i_image_column] == 0 or cloud_mask[i_image_line][i_image_column] == 2 or cloud_mask[i_image_line][i_image_column] == 3) and dust_mask[i_image_line][i_image_column] != 1 );

            std::bitset<Mask::NB_MASK_VALUES> maskValue;
			if(isCloud) 												maskValue.set(Mask::CLOUDS);
			if(dust_mask[i_image_line][i_image_column] == 1) 			maskValue.set(Mask::DUST);
			//cgs 5/17: volcanic ash will be treated in the same way as dust (issue 52) -> changed Mask::VOLCANIC_PLUME to Mask::DUST
			if(volcanic_plume_mask[i_image_line][i_image_column] == 1) 	maskValue.set(Mask::DUST);
			if(isSnow) 													maskValue.set(Mask::SNOW);
			if(isFire)													maskValue.set(Mask::FIRE);
			//if(smoke_mask[i_image_line][i_image_column])				maskValue.set(Mask::SMOKE);

			mask_data.mask(i_line, i_column) = Mask::bitset_to_uchar(maskValue);
		}
	}

}


