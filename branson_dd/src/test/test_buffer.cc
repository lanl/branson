//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_buffer.cc
 * \author Alex Long
 * \date   January 14 2016
 * \brief  Test construction, fill and get functions
 * \note   ***COPYRIGHT_GOES_HERE****
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>

#include "../buffer.h"

int main (void) {

  using std::cout;
  using std::endl;
  using std::vector;

  int nfail = 0;

  // test construction, get and set functions
  {
    bool test_buffer_construction = true;
    Buffer<int> buffer;
    //after construction buffer should have state == EMPTY 
    if (!buffer.empty()) test_buffer_construction = false;
    
    // the vector of objects should also be empty (size == 0)
    if (!buffer.get_object().empty()) test_buffer_construction=false;

    // all other status checks should return false
    if (buffer.received()) test_buffer_construction = false;
    if (buffer.awaiting()) test_buffer_construction = false;
    if (buffer.sent()) test_buffer_construction = false;
    if (buffer.ready()) test_buffer_construction = false;

    if (test_buffer_construction) cout<<"TEST PASSED: Buffer construction"<<endl;
    else { 
      cout<<"TEST FAILED: Buffer construction"<<endl; 
      nfail++;
    }
  }

  // test buffer fill and clear
  {
    bool test_buffer_fill = true;
    Buffer<int> buffer;

    vector<int> test_vector;
    test_vector.push_back(1);
    test_vector.push_back(2);
    test_vector.push_back(3);
    test_vector.push_back(4);

    buffer.fill(test_vector);
 
    //test vector should be inside buffer
    if (buffer.get_object().size() != test_vector.size() )
      test_buffer_fill = false;
    if (buffer.get_object()[0] != test_vector[0]) test_buffer_fill = false;
    if (buffer.get_object()[1] != test_vector[1]) test_buffer_fill = false;
    if (buffer.get_object()[2] != test_vector[2]) test_buffer_fill = false;
    if (buffer.get_object()[3] != test_vector[3]) test_buffer_fill = false;

    //after filling buffer, it should be ready
    if (!buffer.ready()) test_buffer_fill = false;

    // all other status checks should return false
    if (buffer.received()) test_buffer_fill = false;
    if (buffer.awaiting()) test_buffer_fill = false;
    if (buffer.sent()) test_buffer_fill = false;
    if (buffer.empty()) test_buffer_fill = false;

    // clear buffer, buffer should be empty and buffer should have state == EMPTY
    buffer.reset();
    if (!buffer.empty() ) test_buffer_fill = false;
    
    if (test_buffer_fill) cout<<"TEST PASSED: Buffer fill"<<endl;
    else { 
      cout<<"TEST FAILED: Buffer fill"<<endl; 
      nfail++;
    }
  }
  
  return nfail;
}
//---------------------------------------------------------------------------//
// end of test_buffer.cc
//---------------------------------------------------------------------------//
