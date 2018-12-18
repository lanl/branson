//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_counter_rng.cc
 * \author Alex Long
 * \date   April 27, 2017
 * \brief  Simple random number testing
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//
#include <iostream>

#include "../RNG.h"
#include "testing_functions.h"

int main(void) {

  using std::cout;
  using std::endl;

  int nfail = 0;

  // test to make sure rng only produces values in (0,1) and not [0,1] or (0,1] etc.
  {
    bool test_rng_max_min = true;

    double max_rng = 0.0;
    double min_rng = 100.0;

    RNG rng_1;
    rng_1.set_seed(12412);

    double temp_random_double;
    double random_avg = 0.0;
    uint32_t n_sample = 10000000;
    for (uint32_t i = 0; i < n_sample; ++i) {
      temp_random_double = rng_1.generate_random_number();
      random_avg += temp_random_double;
      if (temp_random_double > max_rng)
        max_rng = temp_random_double;
      if (temp_random_double < min_rng)
        min_rng = temp_random_double;
    }
    random_avg = random_avg / n_sample;

    cout.precision(17);
    cout << "Minimum RNG double: " << min_rng << endl;
    cout << "Maximum RNG double: " << max_rng << endl;
    cout << "Average RNG double: " << random_avg << endl;
    cout.precision(6);

    if (max_rng >= 1.0 || max_rng <= 0.0)
      test_rng_max_min = false;
    if (min_rng >= 1.0 || min_rng <= 0.0)
      test_rng_max_min = false;
    if (!soft_equiv(random_avg, 0.5, 1.0e-3))
      test_rng_max_min = false;

    if (test_rng_max_min)
      cout << "TEST PASSED: RNG double in (0,1)" << endl;
    else {
      cout << "TEST FAILED: RNG double not in (0,1)" << endl;
      nfail++;
    }
  }

  // test for identical seed and identical values
  {
    bool test_rng_equal = true;

    RNG rng_1;
    rng_1.set_seed(4056);
    RNG rng_2;
    rng_2.set_seed(4056);

    double temp_random_double_1 = rng_1.generate_random_number();
    double temp_random_double_2 = rng_2.generate_random_number();

    cout << "Double from rng_1: " << temp_random_double_1;
    cout << " Double from rng_2: " << temp_random_double_2 << endl;

    if (temp_random_double_1 != temp_random_double_2)
      test_rng_equal = false;

    temp_random_double_1 = rng_1.generate_random_number();
    temp_random_double_2 = rng_2.generate_random_number();

    cout << "Double from rng_1: " << temp_random_double_1;
    cout << " Double from rng_2: " << temp_random_double_2 << endl;

    if (temp_random_double_1 != temp_random_double_2)
      test_rng_equal = false;

    if (test_rng_equal)
      cout << "TEST PASSED: RNG produces same value with same seed" << endl;
    else {
      cout << "TEST FAILED: RNG does not produce same value with same seed"
           << endl;
      nfail++;
    }
  }

  // test if identical seed and different stream values produces different
  // values
  {
    bool test_rng_stream_equal = true;

    RNG rng_1;
    rng_1.set_seed(4056, 10);
    RNG rng_2;
    rng_2.set_seed(4056, 100);

    double temp_random_double_1 = rng_1.generate_random_number();
    double temp_random_double_2 = rng_2.generate_random_number();

    cout << "Double from rng_1: " << temp_random_double_1;
    cout << " Double from rng_2: " << temp_random_double_2 << endl;

    if (temp_random_double_1 == temp_random_double_2)
      test_rng_stream_equal = false;

    temp_random_double_1 = rng_1.generate_random_number();
    temp_random_double_2 = rng_2.generate_random_number();

    cout << "Double from rng_1: " << temp_random_double_1;
    cout << " Double from rng_2: " << temp_random_double_2 << endl;

    if (temp_random_double_1 == temp_random_double_2)
      test_rng_stream_equal = false;

    if (test_rng_stream_equal) {
      cout << "TEST PASSED: RNG produces different value with same seed and ";
      cout << "different streams" << endl;
    } else {
      cout << "TEST FAILED: RNG produces the same value with the seed and ";
      cout << "different streams" << endl;
      nfail++;
    }
  }

  // test for identical seed and different values
  {
    bool test_rng_diff = true;

    RNG rng_1;
    rng_1.set_seed(1985);
    RNG rng_2;
    rng_2.set_seed(2015);

    double temp_random_double_1 = rng_1.generate_random_number();
    double temp_random_double_2 = rng_2.generate_random_number();

    cout << "Double from rng_1: " << temp_random_double_1;
    cout << " Double from rng_2: " << temp_random_double_2 << endl;

    if (temp_random_double_1 == temp_random_double_2)
      test_rng_diff = false;

    temp_random_double_1 = rng_1.generate_random_number();
    temp_random_double_2 = rng_2.generate_random_number();

    cout << "Double from rng_1: " << temp_random_double_1;
    cout << " Double from rng_2: " << temp_random_double_2 << endl;

    if (temp_random_double_1 == temp_random_double_2)
      test_rng_diff = false;

    if (test_rng_diff) {
      cout << "TEST PASSED: RNG produces different value with different seeds";
      cout << endl;
    } else {
      cout << "TEST FAILED: RNG produces the same value with different seeds";
      cout << endl;
      nfail++;
    }
  }

  return nfail;
}
//----------------------------------------------------------------------------//
// end of test_counter_rng.cc
//----------------------------------------------------------------------------//
