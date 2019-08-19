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
 * AbstractReader.h
 *
 *  Created on: 18 Nov 2015
 *      Author: alix
 */

#ifndef TILE_MAKER_MODULE_READERS_ABSTRACTREADER_H_
#define TILE_MAKER_MODULE_READERS_ABSTRACTREADER_H_

#include "../../Common/Global_header.h"
#include "../../Common/Data.h"
#include "../../Common/Logger.h"
#include "../../Common/TileDimensions.h"


/**
 * \ingroup Tile_Maker
 *
 * \brief   Small abstract class containing a pointer to the Python logger. All Reader classes must inherit from it.
 */
class AbstractReader
{
protected:
	//! Python logger
	std::shared_ptr<AbstractLogger> logger;

	//! flag to check if external files to be read are open
	bool isReady = false;

	AbstractReader(std::shared_ptr<AbstractLogger> logger)
	{
		this->logger = logger;
	}

	virtual ~AbstractReader() { }

public:
	bool is_ready() { return isReady; }

};

#endif /* TILE_MAKER_MODULE_READERS_ABSTRACTREADER_H_ */
