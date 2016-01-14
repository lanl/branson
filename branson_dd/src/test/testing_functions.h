
#include <stdlib.h>

template<typename T>
bool soft_equiv(const T& a, const T& b, double tolerance ) {
  using std::abs;
  return  abs( a-b) < tolerance;
}
