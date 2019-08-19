# See: https://github.com/pybind/cmake_example
#      http://cicero.xyz/v2/remark/github/scisoft/toolbox-talks/master/cmake-kung-fu.mkd

cmake_minimum_required(VERSION 3.0)

project(GEDAP)

#################
# Add languages #
#################

enable_language(CXX)

#########################
# Additional path setup # 
#########################

# Add include path for cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
include(GEDAPUtils)

#########################
# Dependency management #
#########################

# Search for libraries (includes automated Conda-based hinting)
option(ENABLE_AUTO_CONDA_HINTING "Enable Conda-based library path search" ON)

# Deal with fortran specials
if (DEFINED GFORTRAN_LIBRARIES_HINT)
    message(STATUS "Library search hint [libgfortran]: ${GFORTRAN_LIBRARIES_HINT}")
endif ()
find_package(gfortran)

# HDF5
if (ENABLE_AUTO_CONDA_HINTING AND DEFINED ENV{CONDA_PREFIX} AND NOT HDF5_ROOT)
    set(HDF5_ROOT "$ENV{CONDA_PREFIX}")
endif ()
if (DEFINED HDF5_ROOT)
    message(STATUS "Library search hint [HDF5]: ${HDF5_ROOT}")
endif ()
find_package(HDF5 REQUIRED)
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})

# NetCDF
if (ENABLE_AUTO_CONDA_HINTING AND DEFINED ENV{CONDA_PREFIX} AND NOT NETCDF_DIR)
    #message(STATUS "Library search hint [NetCDF] -> using Conda: $ENV{CONDA_PREFIX}")
    set(NETCDF_DIR "$ENV{CONDA_PREFIX}")
endif ()
if (DEFINED NETCDF_DIR)
    message(STATUS "Library search hint [NetCDF]: ${NETCDF_DIR}")
endif ()
find_package(NetCDF REQUIRED)
include_directories(SYSTEM ${NETCDF_INCLUDE_DIRS})

# NetCDF-CXX4
if (ENABLE_AUTO_CONDA_HINTING AND DEFINED ENV{CONDA_PREFIX} AND NOT NetCDF_CXX4_DIR)
    set(NetCDF_CXX4_DIR "$ENV{CONDA_PREFIX}")
endif ()
if (DEFINED NetCDF_CXX4_DIR)
    message(STATUS "Library search hint [NetCDF_CXX4]: ${NetCDF_CXX4_DIR}")
endif ()
find_package(NetCDF_CXX4 REQUIRED)
include_directories(SYSTEM ${NetCDF_CXX4_INCLUDE_DIRS})

# Eigen
if (ENABLE_AUTO_CONDA_HINTING AND DEFINED ENV{CONDA_PREFIX} AND NOT Eigen_DIR)
    #message(STATUS "Library search hint [Eigen] -> using Conda: $ENV{CONDA_PREFIX}/include")
    set(Eigen_DIR "$ENV{CONDA_PREFIX}/include")
endif ()
if (DEFINED Eigen_DIR)
    message(STATUS "Library search hint [Eigen]: ${Eigen_DIR}")
endif ()
find_package(Eigen REQUIRED)
include_directories(SYSTEM ${Eigen_INCLUDE_DIRS})

# pybind11
if (NOT PYTHON_EXECUTABLE)
    if (DEFINED ENV{GEDAP_PYTHON_EXECUTABLE})
        message(STATUS "Library search hint [pybind11] -> using \$ENV{GEDAP_PYTHON_EXECUTABLE}: "
                "$ENV{GEDAP_PYTHON_EXECUTABLE}")
        set(PYTHON_EXECUTABLE "$ENV{GEDAP_PYTHON_EXECUTABLE}")
    elseif (DEFINED ENV{CONDA_PREFIX})
        message(STATUS "Library search hint [pybind11] -> using Conda Python executable: "
                "$ENV{CONDA_PREFIX}/bin/python")
        set(PYTHON_EXECUTABLE "$ENV{CONDA_PREFIX}/bin/python")
    endif ()
endif ()
add_subdirectory(pybind11)
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS} ${PYBIND11_INCLUDE_DIR})

###############
# Build types #
###############

# Usage: cmake -DMyOption=ON MyProjectDirectory
# Options are also displayed when using ccmake and cmake-gui
option(ENABLE_OPENMP "Enable OpenMP support" ON)
option(ENABLE_INSTANTANEOUS "Compile in INSTANTANEOUS mode" ON)
option(ENABLE_TESTING "Compile tests" OFF)

option(ENABLE_SIMPLEPREC "Use simple precision in CISAR" OFF)
mark_as_advanced(ENABLE_SIMPLEPREC)
option(ENABLE_GET_ERROR "Return relative errors for analysis" OFF)
mark_as_advanced(ENABLE_GET_ERROR)

# Declare basic C++ compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -fPIC -fvisibility=hidden")
# message(STATUS "C++ compiler flags: ${CMAKE_CXX_FLAGS}")

# Handle OpenMP option
if (ENABLE_OPENMP)
    message(STATUS "Compiling with OpenMP")
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

# Handle INSTANTANEOUS option
if (ENABLE_INSTANTANEOUS)
    message(STATUS "Compiling in INSTANTANEOUS mode")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DINSTANTANEOUS")
else ()
    message(STATUS "Compiling in TIMESERIES mode")
endif ()

# Handle GET_ERROR option
if (ENABLE_GET_ERROR)
    message(STATUS "Compiling in INSTANTANEOUS mode")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGET_ERROR")
endif ()

######################
# Target definitions #
######################

# Setup output directories
set(INSTALL_DIR ${CMAKE_SOURCE_DIR}/install)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/gedap/cpp)
set(CMAKE_TEST_OUTPUT_DIRECTORY ${INSTALL_DIR}/bin)

# # C++ module: coord_convert
# add_library(coord_convert
#         src/CoordinatesConverter/CoordinatesConverter.cpp)

# C++ module: common
add_library(common INTERFACE)
target_sources(common INTERFACE
        src/Common/Data.h
        src/Commm/Global_header.h
        src/Common/Logger.h
        src/Common/TileDimensions.h)

# C++ module: tile
set(SRC_tile
        src/Tiles/TileNetCDFManager.cpp)
set(H_tile
        src/Tiles/FortranStructures.h
        src/Tiles/TileNetCDFManager.h)

add_library(gedap_tile
        ${SRC_tile} ${H_tile})
target_link_libraries(gedap_tile PUBLIC
        ${NetCDF_CXX4_LIBRARIES})

# C++ module: tile_maker
set(SRC_tile_maker
        src/Tile_Maker_Module/TileMaker.cpp
        src/Tile_Maker_Module/MFGTileMaker.cpp
        src/Tile_Maker_Module/Readers/MaskReader.cpp
        src/Tile_Maker_Module/Readers/MFGMaskReader.cpp
        src/Tile_Maker_Module/Readers/MFGObservationsReader.cpp
        src/Tile_Maker_Module/Readers/ModelParamsReader.cpp
        src/Tile_Maker_Module/Readers/ObservationsReader.cpp
        src/Tile_Maker_Module/Readers/StaticDataReader.cpp)
set(H_tile_maker
        src/Tile_Maker_Module/TileMaker.h
        src/Tile_Maker_Module/MFGTileMaker.h
        src/Tile_Maker_Module/Readers/AbstractReader.h
        src/Tile_Maker_Module/MFGTileMaker.h
        src/Tile_Maker_Module/Readers/MaskReader.h
        src/Tile_Maker_Module/Readers/MFGMaskReader.h
        src/Tile_Maker_Module/Readers/MFGObservationsReader.h
        src/Tile_Maker_Module/Readers/ModelParamsReader.h
        src/Tile_Maker_Module/Readers/ObservationsReader.h
        src/Tile_Maker_Module/Readers/Satellite_Parameters.h
        src/Tile_Maker_Module/Readers/StaticDataReader.h)

add_library(gedap_tile_maker
        ${SRC_tile_maker} ${H_tile_maker})
target_link_libraries(gedap_tile_maker PRIVATE
        ${HDF5_CXX_LIBRARIES}
        gedap_tile)

# C++ module: tile_processor
set(SRC_tile_processor
        src/Tile_Processor_Module/PriorManager.cpp
        src/Tile_Processor_Module/TileProcessor.cpp)
set(H_tile_processor
        src/Tile_Processor_Module/PriorManager.h
        src/Tile_Processor_Module/TileProcessor.h)

add_library(gedap_tile_processor
        ${SRC_tile_processor} ${H_tile_processor})
target_link_libraries(gedap_tile_processor PRIVATE
        gedap_tile
        gedap_tile_maker)

add_subdirectory(src/Fortran)

# Setup GEDAP Python lib compilation
pybind11_add_module(libGEDAP
        src/GEDAP.cpp)
target_link_libraries(libGEDAP PRIVATE
        ${LIBGFORTRAN_LIBRARIES}
        gedap_tile
        gedap_tile_maker
        gedap_tile_processor
        CISAR_gfortran_INSTANTANEOUS)

if (ENABLE_OPENMP)
    target_link_libraries(libGEDAP PRIVATE
            ${OpenMP_Fortran_LIBRARY})
endif()

