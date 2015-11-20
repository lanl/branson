#ifndef request_h_
#define request_h_

#include <boost/mpi.hpp>

namespace mpi = boost::mpi;



class Request {

  public:
  Request() 
  : b_valid(false)
    {}
  ~Request() {}
  
  void request(mpi::request _req) {
    req = _req;
    b_valid = true;
  } 

  bool test(void) {
    if (b_valid) {
      if (req.test()) {
        b_valid = false;
        return true;
      }
    }
    return false;
  }

  void cancel(void) {
    if (b_valid) req.cancel();
  }
 
  bool valid(void) {return b_valid;}

  private:
  bool b_valid;
  mpi::request req;  

};

#endif // request_h_
