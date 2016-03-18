#include <cmath>

#ifndef testing_functions_h_
#define testing_functions_h_

bool soft_equiv(const double& a, const double& b, double tolerance ) {
  double diff = a-b;
  return  std::fabs(diff) < tolerance;
}

#endif // #ifndef testing_functions_h_
