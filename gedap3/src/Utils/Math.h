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
#ifndef GEDAP_UTILS_MATH_H
#define GEDAP_UTILS_MATH_H
namespace GEDAP {
template<typename T>
T sqr(T x) {
    return x*x;
}
template<typename T>
T cubic_interpolation(T p[4], T x[4], T x0) {
    T a = 2.0 * p[1] - 2.0 * p[2] + (p[2] - p[0]) * (x[2] - x[1]) / (x[2] - x[0]) +
          (p[3] - p[1]) * (x[2] - x[1]) / (x[3] - x[1]);
    T b = -3.0 * p[1] + 3.0 * p[2] - 2.0 * (p[2] - p[0]) * (x[2] - x[1]) / (x[2] - x[0]) -
          (p[3] - p[1]) * (x[2] - x[1]) / (x[3] - x[1]);
    T c = (p[2] - p[0]) * (x[2] - x[1]) / (x[2] - x[0]);
    T d = p[1];
    x0 = (x0 - x[1]) / (x[2] - x[1]);
    return d + x0 * (c + x0 * (b + x0 * a));
}
template<typename T>
T bicubic_interpolation(T p[4][4], T x[4][4], T y[4][4], T x0, T y0) {
    T arr[4];
    arr[0] = cubic_interpolation(p[0], y[0], y0);
    arr[1] = cubic_interpolation(p[1], y[1], y0);
    arr[2] = cubic_interpolation(p[2], y[2], y0);
    arr[3] = cubic_interpolation(p[3], y[3], y0);
    return cubic_interpolation(arr, x[0], x0);
}
}
#endif
