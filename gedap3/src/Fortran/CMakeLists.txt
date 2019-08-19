cmake_minimum_required(VERSION 3.0)

project(CISAR_gfortran_INSTANTANEOUS LANGUAGES Fortran)
set(CMAKE_Fortran_FLAGS "-funroll-all-loops -O3 -fPIC -fopenmp")
add_library(CISAR_gfortran_INSTANTANEOUS
        minpack.f90
        library.f90)
