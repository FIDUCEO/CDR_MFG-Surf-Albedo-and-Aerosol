# - FindNetCDF_CXX4
# Find NetCDF_CXX4 includes and library
#
# Search hint may be provided through variable NetCDF_CXX4_DIR.
#
# This module returns these variables for the rest of the project to use.
#
#  NetCDF_CXX4_FOUND        - set to true if NetCDF_CXX4 was found
#  NetCDF_CXX4_INCLUDE_DIRS - the directory containing NetCDF_CXX4 header files
#  NetCDF_CXX4_LIBRARIES    - path to NetCDF_CXX4 library files
#
# Normal usage would be:
#  find_package (NetCDF_CXX4 REQUIRED)
#  target_link_libraries (some_target ${NetCDF_CXX4_LIBRARIES})
#
# Author: Vincent Leroy <vincent.leroy@rayference.eu>
# Inspired by: https://gitlab.com/wxmetvis/met.3d

message(STATUS "Searching for NetCDF_CXX4...")

set(NetCDF_CXX4_CHECK_PATHS 
    "/usr/local" # Homebrew
    "/opt/local" # MacPorts
    "/usr")
    
# Standard stuff did not work, so we try to locate Eigen in include paths
find_path(NetCDF_CXX4_INCLUDE_DIR NAMES netcdf
    HINTS ${NetCDF_CXX4_DIR}
    PATHS ${NetCDF_CXX4_CHECK_PATHS}
    PATH_SUFFIXES include)

if (NetCDF_CXX4_INCLUDE_DIR)
    # We found the include dir, we add it to the standard variable
    # ${<library>_INCLUDE_DIRS}
    set (NetCDF_CXX4_INCLUDE_DIRS ${NetCDF_CXX4_INCLUDE_DIR})
endif (NetCDF_CXX4_INCLUDE_DIR)

# Search for library files
find_library (NetCDF_CXX4_LIBRARY NAMES libnetcdf-cxx4.dylib libnetcdf-cxx4.so libnetcdf-c++4.so
    HINTS "${NetCDF_CXX4_DIR}"
    PATHS "${NetCDF_CXX4_CHECK_PATHS}"
    PATH_SUFFIXES lib)

if (NetCDF_CXX4_LIBRARY)
    set (NetCDF_CXX4_LIBRARIES ${NetCDF_CXX4_LIBRARY})
endif (NetCDF_CXX4_LIBRARY)

# Handle REQUIRED/QUIET optional arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF_CXX4
    REQUIRED_VARS NetCDF_CXX4_LIBRARIES NetCDF_CXX4_INCLUDE_DIRS)
