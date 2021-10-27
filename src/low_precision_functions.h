//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   low_precision_functions.h
 * \author Alex Long
 * \date   August 4 2021
 * \brief  Low precision versions of exp and log
 * \note   Copyright (C) 2021 Triad National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef low_precision_functions_h_
#define low_precision_functions_h_

#include "../singem/singem.h"

#include <math.h>
#include <iostream>

//! Low fp exp is meant for arguments in (0,-9.5), this uses the properties of exponents to multiply
// a known exponent at integer points (-1, -2, ..., -9) by the remainder, i.e. exp(-2.5) is the
// value from a table for exp(-2) multiplied by a taylor series expansion for exp(-0.5)
template <typename T>
inline T low_fp_exp_4(T x) {
  constexpr T  exp_at_integers[] = {1.0,
                                        0.36787944117144233,
                                        0.1353352832366127,
                                        0.049787068367863944,
                                        0.01831563888873418,
                                        0.006737946999085467,
                                        0.0024787521766663585,
                                        0.0009118819655545162,
                                        0.00033546262790251185,
                                        0.00012340980408667956};

  int int_part = abs(static_cast<int>(round(x)));
  int_part = (int_part <= 9) ? int_part : 9;
  T exp_x = exp_at_integers[int_part];
  T remainder = x - round(x);

  // 4th order
  T exp_remainder  = 1.0 +  remainder + 0.5 * remainder*remainder + remainder*remainder*remainder * 0.1666666666666667 + remainder*remainder*remainder*remainder*(0.04166666666666);

  //std::cout<<"x: "<<x<<", exp_x: "<<exp_x<<" remainder: "<<remainder<<" exp_remainder: "<<exp_remainder<<" at integer: "<<static_cast<int>(round(x))<<std::endl;
  return (x < -9.0) ? 0.0 : exp_remainder*exp_x;
}

//! Low fp exp is meant for arguments in (0,-9.5), this uses the properties of exponents to multiply
// a known exponent at integer points (-1, -2, ..., -9) by the remainder, i.e. exp(-2.5) is the
// value from a table for exp(-2) multiplied by a taylor series expansion for exp(-0.5)
template <typename T>
inline T low_fp_exp_3(T x) {
  constexpr T  exp_at_integers[] = {1.0,
                                        0.36787944117144233,
                                        0.1353352832366127,
                                        0.049787068367863944,
                                        0.01831563888873418,
                                        0.006737946999085467,
                                        0.0024787521766663585,
                                        0.0009118819655545162,
                                        0.00033546262790251185,
                                        0.00012340980408667956};

  int int_part = abs(static_cast<int>(round(x)));
  int_part = (int_part <= 9) ? int_part : 9;
  T exp_x = exp_at_integers[int_part];
  T remainder = x - round(x);

  // 3rd order
  T exp_remainder  = 1.0 +  remainder + 0.5 * remainder*remainder + remainder*remainder*remainder * 0.1666666666666667;

  //std::cout<<"x: "<<x<<", exp_x: "<<exp_x<<" remainder: "<<remainder<<" exp_remainder: "<<exp_remainder<<" at integer: "<<static_cast<int>(round(x))<<std::endl;
  return (x < -9.0) ? 0.0 : exp_remainder*exp_x;
}

//! Low fp exp is meant for arguments in (0,-9.5), this uses the properties of exponents to multiply
// a known exponent at integer points (-1, -2, ..., -9) by the remainder, i.e. exp(-2.5) is the
// value from a table for exp(-2) multiplied by a taylor series expansion for exp(-0.5)
template <typename T>
inline T simple_low_fp_exp_3(T x) {

  T closest_int;
  T int_part;

  if (x > -0.5) {
    closest_int = 0.0;
    int_part = 1.0;
  }
  else if (x > -1.5) {
    closest_int = 1.0;
    int_part = 0.36787944117144233;
  }
  else if (x > -2.5) {
    closest_int = 2.0;
    int_part = 0.1353352832366127;
  }
  else if (x > -3.5) {
    closest_int = 3.0;
    int_part = 0.049787068367863944;
  }
  else if (x > -4.5) {
    closest_int = 4.0;
    int_part =  0.01831563888873418;
  }
  else if (x > -5.5) {
    closest_int = 5.0;
    int_part = 0.006737946999085467;
  }
  else if (x > -6.5) {
    closest_int = 6.0;
    int_part = 0.0024787521766663585;
  }
  else if (x > -7.5) {
    closest_int = 7.0;
    int_part = 0.0009118819655545162;
  }
  else {
    closest_int = 8.0;
    int_part = 0.00033546262790251185;
  }

  T remainder = x + closest_int;

  // 3rd order
  T exp_remainder  = 1.0 +  remainder + 0.5 * remainder*remainder + remainder*remainder*remainder * 0.1666666666666667;

  return (x < -9.0) ? 0.0 : exp_remainder*int_part;
}


//! A version that converts double to scFloat, does the math and converts back
inline double double_approx_low_fp_exp_4(double x_double) {

  scFloat x(x_double);
  scFloat closest_int;
  scFloat int_part;

  if (x > scFloat(-0.5)) {
    closest_int = 0.0;
    int_part = 1.0;
  }
  else if (x > scFloat(-1.5)) {
    closest_int = 1.0;
    int_part = 0.36787944117144233;
  }
  else if (x > scFloat(-2.5)) {
    closest_int = 2.0;
    int_part = 0.1353352832366127;
  }
  else if (x > scFloat(-3.5)) {
    closest_int = 3.0;
    int_part = 0.049787068367863944;
  }
  else if (x > scFloat(-4.5)) {
    closest_int = 4.0;
    int_part =  0.01831563888873418;
  }
  else if (x > scFloat(-5.5)) {
    closest_int = 5.0;
    int_part = 0.006737946999085467;
  }
  else if (x > scFloat(-6.5)) {
    closest_int = 6.0;
    int_part = 0.0024787521766663585;
  }
  else if (x > scFloat(-7.5)) {
    closest_int = 7.0;
    int_part = 0.0009118819655545162;
  }
  else {
    closest_int = 8.0;
    int_part = 0.00033546262790251185;
  }

  scFloat remainder = x + closest_int;

  // 3rd order
  //scFloat exp_remainder  = scFloat(1.0) +  remainder + scFloat(0.5) * remainder*remainder + remainder*remainder*remainder * scFloat(0.1666666666666667);
  // 4th  order
  scFloat exp_remainder  = scFloat(1.0) +  remainder + scFloat(0.5) * remainder*remainder + remainder*remainder*remainder * scFloat(0.1666666666666667) + remainder*remainder*remainder*remainder*scFloat(0.04166666666666);

  scFloat result = (x < scFloat(-9.0)) ? scFloat(0.0) : exp_remainder*int_part;
  return result.as_double();
}


//! Low fp exp is meant for arguments in (0,-9.5), this uses the properties of exponents to multiply
// a known exponent at integer points (-1, -2, ..., -9) by the remainder, i.e. exp(-2.5) is the
// value from a table for exp(-2) multiplied by a taylor series expansion for exp(-0.5)
inline scFloat approx_low_fp_exp_3(scFloat x) {

  scFloat closest_int;
  scFloat int_part;

  if (x > scFloat(-0.5)) {
    closest_int = 0.0;
    int_part = 1.0;
  }
  else if (x > scFloat(-1.5)) {
    closest_int = 1.0;
    int_part = 0.36787944117144233;
  }
  else if (x > scFloat(-2.5)) {
    closest_int = 2.0;
    int_part = 0.1353352832366127;
  }
  else if (x > scFloat(-3.5)) {
    closest_int = 3.0;
    int_part = 0.049787068367863944;
  }
  else if (x > scFloat(-4.5)) {
    closest_int = 4.0;
    int_part =  0.01831563888873418;
  }
  else if (x > scFloat(-5.5)) {
    closest_int = 5.0;
    int_part = 0.006737946999085467;
  }
  else if (x > scFloat(-6.5)) {
    closest_int = 6.0;
    int_part = 0.0024787521766663585;
  }
  else if (x > scFloat(-7.5)) {
    closest_int = 7.0;
    int_part = 0.0009118819655545162;
  }
  else {
    closest_int = 8.0;
    int_part = 0.00033546262790251185;
  }

  scFloat remainder = x + closest_int;

  // 3rd order
  //scFloat exp_remainder  = scFloat(1.0) +  remainder + scFloat(0.5) * remainder*remainder + remainder*remainder*remainder * scFloat(0.1666666666666667);
  // 4th  order
  scFloat exp_remainder  = scFloat(1.0) +  remainder + scFloat(0.5) * remainder*remainder + remainder*remainder*remainder * scFloat(0.1666666666666667) + remainder*remainder*remainder*remainder*scFloat(0.04166666666666);

  return (x < scFloat(-9.0)) ? scFloat(0.0) : exp_remainder*int_part;
}

template <typename T>
T log_taylor(T x) {
  value = 0.0
  // the mercater series gives the solution for log(1 - x), so if you want the the log of x, do
  // taylor_series_log(1-x)
  T x_trans = 1.0 - x;
  return x_trans + 0.5*x_trans*x_trans + 0.3333333333*x_trans*x_trans
  /*
  for i in range(1, order+1):
    value += x_trans**(i)/i
  return -value
  */
}

template <typename T>
T log_taylor_div_actual(T x) {
 T  value = 0.0
  if (x < 0.03125)
    return log_taylor(x*32) - log(32)
  else if (x < 0.0625)
    return log_taylor(x*16) - log(16)
  else if (x < 0.125)
    return log_taylor(x*8) - log(8)
  else if (x < 0.25)
    return log_taylor(x*4) - log(4)
  else if(x< 0.5)
    return log_taylor(x*2) - log(2)
  else if(x< 0.66666666):
    return log_taylor(x*1.5) - log(1.5)
  else
    return log_taylor(x)
}

//! Low fp log is meant for arguments in (0, 1), this uses the properties of log to multiply
// a known exponent at integer points (-1, -2, ..., -9) by the remainder, i.e. exp(-2.25) is the
// value from a table for exp(-2) multiplied by a taylor series expansion for exp(-0.25)
/*
inline double low_fp_exp(double x) {
  constexpr double exp_at_integers[] = {1.0,
                                        0.36787944117144233,
                                        0.1353352832366127,
                                        0.049787068367863944,
                                        0.01831563888873418,
                                        0.006737946999085467,
                                        0.0024787521766663585,
                                        0.0009118819655545162,
                                        0.00033546262790251185,
                                        0.00012340980408667956};

  double exp_x = exp_at_integers[abs(static_cast<int>(round(x)))];
  double remainder = x - round(x);
  double exp_remainder  = 1.0 +  remainder + 0.5 * remainder*remainder + remainder*remainder * remainder * 0.1666666666666667;
  //std::cout<<"x: "<<x<<", exp_x: "<<exp_x<<" remainder: "<<remainder<<" exp_remainder: "<<exp_remainder<<" at integer: "<<static_cast<int>(round(x))<<std::endl;
  return (x < -9.3) ? 0.0 : exp_remainder*exp_x;
}
*/


#endif
//---------------------------------------------------------------------------//
// end of low_precision_functions.h
//---------------------------------------------------------------------------//

