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
 * CloudMaskReader.h
 *
 *  Created on: 13 Oct 2015
 *      Author: alix
 */

#ifndef TILE_MAKER_MODULE_CLOUDMASKREADER_H_
#define TILE_MAKER_MODULE_CLOUDMASKREADER_H_

#include <H5Cpp.h>

#include "AbstractReader.h"
#include "Satellite_Parameters.h"

//#define DEBUG_MASK


/**
 * \ingroup Tile_Maker
 *
 * \brief   Read external masks and aggregate them in one unique mask. The data are saved in an instance of Mask structure.
 * \details This class use the STL bitset class. \sa Mask
 *
 * \warning Fire and Jon Snow are currently not taken into account.
 * \warning External data are assumed to be already in the right projection.
 *
 * \todo    Modify the class to include the fire and snow.
 */
class MaskReader : public AbstractReader
{
	std::string inputFile;
	std::shared_ptr<H5::H5File> HDFinputFile = nullptr;

	unsigned char  cloud_mask[IMAGE_NB_LINES][IMAGE_NB_COLUMNS];
	unsigned char  dust_mask[IMAGE_NB_LINES][IMAGE_NB_COLUMNS];
	unsigned char  volcanic_plume_mask[IMAGE_NB_LINES][IMAGE_NB_COLUMNS];
	unsigned char  fire_mask[IMAGE_NB_LINES][IMAGE_NB_COLUMNS];
	bool           smoke_mask[IMAGE_NB_LINES][IMAGE_NB_COLUMNS];
	unsigned short test_mask[IMAGE_NB_LINES][IMAGE_NB_COLUMNS];

public:
	MaskReader(const std::string inputFile, std::shared_ptr<AbstractLogger> logger);

	virtual ~MaskReader();

	/**
	 * Reads external mask files and store the data
	 */
	void read();

	/**
	 * For a given tile, aggregates the external masks into one and copy it into an instance of Mask structure
	 */
	void extract_data(Mask& mask_data, TileDimensions& tile_dims);

private:
	//! For debugging
	void print_buffer(std::string name, unsigned char** a_buffer, const int nbLines, const int nbColumns);
	//! For debugging
	void print_buffer(std::string name, unsigned short** a_buffer, const int nbLines, const int nbColumns);
};

#endif /* TILE_MAKER_MODULE_CLOUDMASKREADER_H_ */
