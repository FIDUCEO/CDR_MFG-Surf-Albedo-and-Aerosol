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
 * ObservationsReader.h
 *
 *  Created on: 13 Oct 2015
 *      Author: alix
 */

#ifndef TILE_MAKER_MODULE_OBSERVATIONSREADER_H_
#define TILE_MAKER_MODULE_OBSERVATIONSREADER_H_

#include <map>
#include "AbstractReader.h"


/**
 * \ingroup Tile_Maker
 *
 * \brief   Abstract class definning 2 pure virtual methods: read and extract_data.
 * \details Any class designed to read observations data MUST inherit from this class and implement the 2 virtual methods.
 *
 * \todo    For each new radiometer or satellite to include in GEDAP, create a new class that inherit from this one and
 *          implement the virtual method extract_data.
 *          See for example the class MSG2ObservationsReader.
 */
class ObservationsReader : public AbstractReader {

public:
    explicit ObservationsReader(std::shared_ptr<AbstractLogger> logger)
        : AbstractReader(logger) {}

    virtual ~ObservationsReader() = default;

    virtual void read() = 0;

    virtual void extract_data(Observations &obs_data, TileDimensions &tile_dims) = 0;
};

#endif /* TILE_MAKER_MODULE_OBSERVATIONSREADER_H_ */
