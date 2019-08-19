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
#ifndef SOURCE_CPP_TILE_MAKER_MODULE_TILEMAKER_GLOBAL_HEADER_H_
#define SOURCE_CPP_TILE_MAKER_MODULE_TILEMAKER_GLOBAL_HEADER_H_
#define SUCCESS 0
#define ERROR 1
#define ERROR_OPEN_HDF5_FILE 2
#define ERROR_READ_HDF5_FILE 3
#define ERROR_WRIE_HDF5_FILE 4
#define ERROR_OPEN_NETCDF_FILE 5
#define ERROR_READ_NETCDF_FILE 6
#define ERROR_WRIE_NETCDF_FILE 7
#define PI (double)(3.14159265358979323846)
#define TWO_PI (double)(2.0 * PI)
#define EPSILON ((double)1.0E-10)
#define MAX_NB_LINES 400
#define MAX_NB_COLUMNS 400
#define MAX_NB_PIXELS (MAX_NB_LINES * MAX_NB_COLUMNS)
#define MAX_NB_BANDS 1
#define MAX_NB_AEROSOLS_CLASSES 1
#define MAX_NB_ANGLES 64
#define QUADRATURE_POINTS 8
#define MAX_NB_ATM 1
#define NB_SURFACE_PARAMS 4
#define MAX_NB_OBS (MAX_NB_BANDS * MAX_NB_ANGLES)
#define MAX_NB_AOT_VALUES (MAX_NB_BANDS * MAX_NB_AEROSOLS_CLASSES * MAX_NB_ANGLES)
#define MAX_NB_VAR_PER_BAND (NB_SURFACE_PARAMS + (MAX_NB_ANGLES * MAX_NB_AEROSOLS_CLASSES))
#define MAX_NB_VAR (MAX_NB_BANDS * MAX_NB_VAR_PER_BAND)
#define MAX_ARRAY_SIZE MAX_NB_VAR + MAX_NB_OBS
#define MAX_NB_STATE_VECTOR (NB_SURFACE_PARAMS + MAX_NB_AEROSOLS_CLASSES)
#define NB_QUALITY_TEST 7
#define COEFF_DEG2RAD (double)(PI/180.0)
#define COEFF_RAD2DEG (double)(180./PI)
#define COEFF_DU2KGM2 (double)(2.1415e-5)
#define COEFF_KGM22DU (1./COEFF_DU2KGM2)
enum RPV_Params_Order {
    RHO_0, K, THETA, RHO_C
};
constexpr static double MISSING_VALUE = -999.0;
constexpr static double MIN_SZA = 0.0;
constexpr static double MAX_SZA = 70.0;
constexpr static double MIN_VZA = 0.0;
constexpr static double MAX_VZA = 70.0;
constexpr static double SZA_THRESHOLD = 70.0;
constexpr static double MIN_RHO_0 = 1.e-8;
constexpr static double MIN_K = 1.e-8;
constexpr static double MIN_THETA = -0.5 + 1.e-8;
constexpr static double MIN_RHO_C = 1.e-8;
constexpr static double MIN_AOT = 1.e-8;
constexpr static double MAX_RHO_0 = 1.0;
constexpr static double MAX_K = 1.0;
constexpr static double MAX_THETA = 1.0;
constexpr static double MAX_RHO_C = 1.0;
constexpr static double MAX_AOT = 5.0;
constexpr static double MIN_RPV_VALUES[NB_SURFACE_PARAMS] = {MIN_RHO_0, MIN_K, MIN_THETA, MIN_RHO_C};
constexpr static double MAX_RPV_VALUES[NB_SURFACE_PARAMS] = {MAX_RHO_0, MAX_K, MAX_THETA, MAX_RHO_C};
#include <ctime>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <bitset>
#include <thread>
#include <memory>
#include <map>
#include <assert.h>
#include <math.h>
#include <eigen3/Eigen/Dense>
namespace Eigen {
    typedef Matrix<unsigned char, Dynamic, Dynamic, RowMajor> MatrixXuc_RowMaj;
    typedef Matrix<int, Dynamic, Dynamic, RowMajor> MatrixXi_RowMaj;
    typedef Matrix<float, Dynamic, Dynamic, RowMajor> MatrixXf_RowMaj;
    typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXd_RowMaj;
    typedef Array<int, Dynamic, RowMajor> ArrayXi_RowMaj;
    typedef Array<float, Dynamic, RowMajor> ArrayXf_RowMaj;
    typedef Array<double, Dynamic, RowMajor> ArrayXd_RowMaj;
    typedef Array<double, 4, 4, RowMajor> Array44d_RowMaj;
    typedef Array<short, Dynamic, Dynamic, RowMajor> ArrayXXs_RowMaj;
    typedef Array<unsigned short, Dynamic, Dynamic, RowMajor> ArrayXXus_RowMaj;
    typedef Array<unsigned char, Dynamic, Dynamic, RowMajor> ArrayXXuc_RowMaj;
    typedef Array<signed char, Dynamic, Dynamic, RowMajor> ArrayXXsc_RowMaj;
    typedef Array<int, Dynamic, Dynamic, RowMajor> ArrayXXi_RowMaj;
    typedef Array<unsigned int, Dynamic, Dynamic, RowMajor> ArrayXXui_RowMaj;
    typedef Array<float, Dynamic, Dynamic, RowMajor> ArrayXXf_RowMaj;
    typedef Array<double, Dynamic, Dynamic, RowMajor> ArrayXXd_RowMaj;
    typedef Array<float, Dynamic, Dynamic, ColMajor> ArrayXXf_ColMaj;
}
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include "../Utils/Time.h"
inline double deg_to_rad(double deg) { return deg * COEFF_DEG2RAD; }
inline double rad_to_deg(double rad) { return rad * COEFF_RAD2DEG; }
inline void check_RPV_physical_range(int param_index, double &value) {
    if (value < MIN_RPV_VALUES[param_index]) value = MIN_RPV_VALUES[param_index];
    if (value > MAX_RPV_VALUES[param_index]) value = MAX_RPV_VALUES[param_index];
}
inline void check_AOT_physical_range(double &value) {
    if (value < MIN_AOT) value = MIN_AOT;
    if (value > MAX_AOT) value = MAX_AOT;
}
inline std::string channel_to_string(const int channel) {
    return (channel < 10) ? "0" + std::to_string(channel) : std::to_string(channel);
}
enum PROCESSING_ERROR {
    None, UnderIlluminated, OutOfView, SunGlint, SunGlintAffected, NegativeBRF, CloudPresence, BadPixelQuality
};
struct Error {
    int nb_lines;
    int nb_columns;
    Eigen::MatrixXf_RowMaj time;
    PROCESSING_ERROR **error;
    PROCESSING_ERROR *storage;
    Error() : Error(0, 0) {}
    Error(const int nb_lines, const int nb_columns) : nb_lines(nb_lines), nb_columns(nb_columns) {
        time = Eigen::MatrixXf_RowMaj::Constant(nb_lines, nb_columns, MISSING_VALUE);
        error = new PROCESSING_ERROR *[nb_lines];
        storage = new PROCESSING_ERROR[nb_lines * nb_columns];
        for (int i = 0; i < nb_lines; ++i) {
            error[i] = storage + nb_columns * i;
        }
    }
    ~Error() {}
    void reset() {
        time.setConstant(MISSING_VALUE);
    }
};
enum class AEROSOL_TYPE {
    FINE_MODE, COARSE_MODE
};
struct Aerosol_Class {
    std::string name;
    AEROSOL_TYPE type;
    Aerosol_Class(std::string name, AEROSOL_TYPE type) : name(name), type(type) {};
    virtual ~Aerosol_Class() {};
};
inline int get_nb_coarse_aer_classes(const std::vector<Aerosol_Class> &v_aerosol_classes) {
    int nb_coarse_aer_classes = 0;
    for (const Aerosol_Class &aer_class : v_aerosol_classes)
        if (aer_class.type == AEROSOL_TYPE::COARSE_MODE)nb_coarse_aer_classes++;
    return nb_coarse_aer_classes;
}
#endif /* SOURCE_CPP_TILE_MAKER_MODULE_TILEMAKER_GLOBAL_HEADER_H_ */
