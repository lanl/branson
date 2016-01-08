#ifndef buffer_h_
#define buffer_h_

#include <vector>

namespace mpi = boost::mpi;


template<class T>
class Buffer {

  public:
  Buffer() 
  : status(EMPTY)
    {}
  ~Buffer() {}
 
  void fill (std::vector<T> _object) {
    object = _object;
    status = READY;
  } 

  std::vector<T>& get_buffer (void) {return object;}

  void reset(void) {object.clear(); status=EMPTY;} 
  void set_sent(void) {status = SENT;}
  void set_received(void) {status = RECEIVED;}
  void set_awaiting(void) {status = AWAITING;}

  bool sent(void) {return status == SENT;}
  bool awaiting(void) {return status == AWAITING;}
  bool ready(void) {return status == READY;}
  bool received(void) {return status == RECEIVED;}
  bool empty(void) {return status == EMPTY;}


  private:
  std::vector<T> object;
  unsigned int status;
  enum {EMPTY, READY, SENT, AWAITING, RECEIVED};

};

#endif // buffer_h_
