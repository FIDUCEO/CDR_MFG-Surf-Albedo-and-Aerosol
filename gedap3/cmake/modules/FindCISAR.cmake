# - FindCISAR
# Find CISAR Fortran library
#
# Search hint **must** be provided through variable CISAR_DIR.
#
# This module returns these variables for the rest of the project to use.
#
#  CISAR_FOUND                   - set to true if CISAR headers and at least one lib were found
#  CISAR_FOUND_LIB               - set to true if at least one lib were found
#  CISAR_FOUND_INSTANTANEOUS     - set to true if CISAR INSTANTANEOUS was found
#  CISAR_FOUND_TIMESERIES        - set to true if CISAR TIMESERIES was found
#  CISAR_INCLUDE_DIRS            - the directory containing the CISAR header files
#  CISAR_LIBRARIES_INSTANTANEOUS - path to CISAR INSTANTANEOUS libraries
#  CISAR_LIBRARIES_TIMESERIES    - path to CISAR TIMESERIES libraries
#
# Author: Vincent Leroy <vincent.leroy@rayference.eu>

message(STATUS "Searching for CISAR...")

# set(CISAR_CHECK_PATHS
#   # Standard CISAR install paths
#   )

# Standard can't not work, so we try to locate CISAR in include paths
find_path(CISAR_INCLUDE_DIR NAMES config.in 
    HINTS ${CISAR_DIR}
    PATHS ${CISAR_CHECK_PATHS}
    PATH_SUFFIXES CONFIG)

if (CISAR_INCLUDE_DIR)
    # We found the include dir, we add it to the standard variable
    # ${<library>_INCLUDE_DIRS}
    set (CISAR_INCLUDE_DIRS ${CISAR_INCLUDE_DIR})
endif (CISAR_INCLUDE_DIR)

# Search for library files
find_library(CISAR_LIBRARIES_INSTANTANEOUS NAMES libCISAR_gfortran_INSTANTANEOUS.a
    HINTS ${CISAR_DIR}
    PATHS ${CISAR_CHECK_PATHS}
    PATH_SUFFIXES LIB)

if (CISAR_LIBRARIES_INSTANTANEOUS)
    set(CISAR_FOUND_INSTANTANEOUS TRUE)
else ()
    message(STATUS "Warning: CISAR INSTANTANEOUS library not found")
endif ()

find_library(CISAR_LIBRARIES_TIMESERIES NAMES libCISAR_gfortran_TIMESERIES.a
        HINTS ${CISAR_DIR}
        PATHS ${CISAR_CHECK_PATHS}
        PATH_SUFFIXES LIB)

if (CISAR_LIBRARIES_TIMESERIES)
    set(CISAR_FOUND_TIMESERIES TRUE)
else ()
    message(STATUS "Warning: CISAR TIMESERIES library not found")
endif ()

if (CISAR_LIBRARIES_INSTANTANEOUS OR CISAR_LIBRARIES_TIMESERIES)
    set(CISAR_FOUND_LIB TRUE)
endif ()

# Handle REQUIRED/QUIET optional arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CISAR
        REQUIRED_VARS CISAR_INCLUDE_DIRS CISAR_FOUND_LIB)
