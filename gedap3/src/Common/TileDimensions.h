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
#ifndef TILEPROPERTIES_H_
#define TILEPROPERTIES_H_
#include <H5Cpp.h>
#include <iostream>
typedef std::tuple<std::string,int,int,int,int> TileDimensions_Tuple;
struct TileDimensions
{
 std::string st_ID = "";
 int first_line = 0;
 int last_line = 0;
 int first_column = 0;
 int last_column = 0;
 int nb_lines = 0;
 int nb_columns = 0;
 int nb_pixels = 0;
public:
 TileDimensions(){ }
 TileDimensions(const std::string st_ID, const int first_line, const int last_line, const int first_column, const int last_column) :
  st_ID(st_ID), first_line(first_line), last_line(last_line), first_column(first_column), last_column(last_column)
 {
  nb_lines = get_nb_lines();
  nb_columns = get_nb_columns();
  nb_pixels = nb_lines * nb_columns;
 }
 int get_nb_lines() { return last_line - first_line + 1; }
 int get_nb_columns() { return last_column - first_column + 1; }
 int get_nb_pixels() { return get_nb_lines() * get_nb_columns(); }
 int get_pixel_index(const int line, const int column) { return line * nb_lines + column; }
 void set_start_count_vectors(std::vector<size_t>& start, std::vector<size_t>& count)
 {
  start.clear();
  start.push_back(first_line);
  start.push_back(first_column);
  count.push_back(nb_lines);
  count.push_back(nb_columns);
 }
 void set_offset_dims_arrays(hsize_t offset[2], hsize_t dims[2])
 {
  offset[0] = first_line;
  offset[1] = first_column;
  dims[0] = nb_lines;
  dims[1] = nb_columns;
 }
    void dump()
    {
  std::cout << "TILE PROPERTIES: " << std::endl;
  std::cout << "index: " << st_ID << std::endl;
  std::cout << "first line, last line, first column, last column: " << first_line << " " << last_line << " "
      << first_column << " " << last_column << std::endl;
  std::cout << "nb lines, columns, pixels: " << nb_lines << " " << nb_columns << " " << nb_pixels << std::endl;
  std::cout << std::endl;
    }
};
#endif
