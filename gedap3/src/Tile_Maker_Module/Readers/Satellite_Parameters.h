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
 * Satellite_Parameters.h
 *
 *  Created on: 11 Dec 2015
 *      Author: alix
 */

#ifndef TILE_MAKER_MODULE_READERS_MSG2_PARAMETERS_H_
#define TILE_MAKER_MODULE_READERS_MSG2_PARAMETERS_H_

/**
 * \ingroup Tile_Maker
 * \file    Satellite_Parameters.h
 */

#define IMAGE_NB_BANDS      1

// Size for SEVIRI
//#define IMAGE_NB_LINES      3712
//#define IMAGE_NB_COLUMNS    3712

// Size for MVIRI
#define IMAGE_NB_LINES      5000
#define IMAGE_NB_COLUMNS    5000
#define ERA_INTERIM_IMAGE_NB_LINES    241
#define ERA_INTERIM_IMAGE_NB_COLUMNS  480

#define IMAGE_NB_PIXELS     (IMAGE_NB_LINES * IMAGE_NB_COLUMNS)
#define MSG2_NB_BITS_COUNTS 10

#endif /* TILE_MAKER_MODULE_READERS_MSG2_PARAMETERS_H_ */
