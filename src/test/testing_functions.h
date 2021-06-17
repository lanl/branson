//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   testing_functions.h
 * \author Alex Long
 * \date   December 10 2015
 * \brief  Provide soft equivalence functions for testing
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//
#include <cmath>

#ifndef testing_functions_h_
#define testing_functions_h_

bool soft_equiv(const double& a, const double& b, double tolerance ) {
  double diff = a-b;
  return  std::fabs(diff) < tolerance;
}

#endif // testing_functions_h_
//---------------------------------------------------------------------------//
// end of testing_functions.h
//---------------------------------------------------------------------------//

