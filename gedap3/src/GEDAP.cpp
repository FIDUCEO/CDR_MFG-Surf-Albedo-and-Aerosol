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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
#include "Common/TileDimensions.h"
#include "Tile_Maker_Module/MFGTileMaker.h"
#include "Tile_Processor_Module/TileProcessor.h"
#include "Common/Logger.h"
class PythonMFGTileMaker : public MFGTileMaker {
public:
 PythonMFGTileMaker(const std::vector<std::string> &v_bands,
        const std::vector<TileDimensions_Tuple> &v_tiles_as_tuples,
        py::object py_logger)
   : MFGTileMaker(
   v_bands,
   v_tiles_as_tuples,
   std::make_shared<PythonLogger>(py_logger)) {
 }
 virtual ~PythonMFGTileMaker() {}
};
class PythonTileProcessor : public TileProcessor {
public:
    PythonTileProcessor(
   const std::string cisar_lut_directory, const std::string cisar_lut_ID,
   const std::string cisar_aeromod_filename, const int max_nb_iterations, const int nb_days_acc_period,
   const std::vector<std::string> &v_bands, const std::vector<std::string> &map_aerosol_classes,
   const bool b_temporal_smoothness, const bool b_spectral_constraint,
   const bool b_surface_spectral_constraint, py::object py_logger)
   : TileProcessor(
   cisar_lut_directory,
   cisar_lut_ID,
   cisar_aeromod_filename,
   max_nb_iterations,
   nb_days_acc_period,
   v_bands,
   map_aerosol_classes,
   b_temporal_smoothness,
   b_spectral_constraint,
   b_surface_spectral_constraint,
   std::make_shared<PythonLogger>(py_logger)) {
 }
 virtual ~PythonTileProcessor() { }
};
PYBIND11_MODULE(libGEDAP, m) {
    m.doc() = "C++ module for the creation of static and input tiles";
    py::class_<PythonMFGTileMaker>(m, "MFGTileMaker")
            .def(py::init<const std::vector<std::string>&, const std::vector<TileDimensions_Tuple>&, py::object>())
            .def("make_static_data_tile", &MFGTileMaker::make_static_data_tile)
            .def("make_input_tile", &MFGTileMaker::make_input_tile)
            .def("read_first_ERA_INTERIM", &MFGTileMaker::read_first_ERA_INTERIM)
   ;
    py::class_<PythonTileProcessor>(m, "TileProcessor")
   .def(py::init<const std::string, const std::string, const std::string, const int, const int,
           const std::vector<std::string> &, const std::vector<std::string> &,
           const bool, const bool, const bool, py::object>())
   .def("prepare_to_process_new_tile", &TileProcessor::prepare_to_process_new_tile)
   .def("set_params_sigmoid", &TileProcessor::set_params_sigmoid)
   .def("set_params_log", &TileProcessor::set_params_log)
   .def("set_sigma_spectral_constraint", &TileProcessor::set_sigma_spectral_constraint)
   .def("load_static_data_tile", &TileProcessor::load_static_data_tile)
   .def("load_input_tiles", &TileProcessor::load_input_tiles)
   .def("load_prior_tile", &TileProcessor::load_prior_tile)
   .def("create_default_prior_tile", &TileProcessor::create_default_prior_tile)
   .def("inversion", &TileProcessor::inversion)
   ;
 py::register_exception<CisarError>(m, "CisarError");
}
