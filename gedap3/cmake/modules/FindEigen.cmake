# - FindEigen
# Find Eigen header-only library
#
# Search hint may be provided through variable NetCDF_CXX4_DIR.
#
# This module sets the following variables:
#
#  Eigen_FOUND        - set to true if Eigen was found
#  Eigen_INCLUDE_DIRS - the directory containing the Eigen header files
#
# To give an hint for the directory to look for the user must set Eigen_DIR
#
# Author: Vincent Leroy <vincent.leroy@rayference.eu>
# Original author: Ruben Di Battista <rubendibattista@gmail.com>

message(STATUS "Searching for Eigen...")

set(Eigen_CHECK_PATHS 
    "/usr/local/include"
    "/usr/local/homebrew/include" # Mac OS X
    "/opt/local/var/macports/software" # Mac OS X.
    "/opt/local/include"
    "/usr/include")

# We try to locate Eigen in include paths
find_path(Eigen_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
    HINTS ${Eigen_DIR}
    PATHS ${Eigen_CHECK_PATHS}
    PATH_SUFFIXES eigen3 Eigen3 eigen Eigen)

if (Eigen_INCLUDE_DIR)
    # We found the include dir, we add it to the standard variable
    # ${<library>_INCLUDE_DIRS}
    set (Eigen_INCLUDE_DIRS ${Eigen_INCLUDE_DIR})
endif (Eigen_INCLUDE_DIR)

# Handle REQUIRED/QUIET optional arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen
  REQUIRED_VARS Eigen_INCLUDE_DIRS)
