/*
  Author: Alex Long
  Date: 12/15/2015
  Name: buffer.h
*/
#ifndef buffer_h_
#define buffer_h_

#include <vector>

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

  void *get_buffer (void) {
    if (object.size()>0) return &object[0];
    else return &object;
  }
  
  std::vector<T>& get_object(void) {return object;}
  
  void resize(uint32_t new_size) {object.resize(new_size);}
  void clear(void) {object.clear(); status=EMPTY;} 
  void reset(void) {status=EMPTY;} 
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
  uint32_t status;
  enum {EMPTY, READY, SENT, AWAITING, RECEIVED};

};

#endif // buffer_h_
