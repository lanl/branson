//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   buffer.h
 * \author Alex Long
 * \date   December 15 2015
 * \brief  Class for simplifying data management in MPI messaging
 * \note   ***COPYRIGHT_GOES_HERE****
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

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

  bool sent(void) const {return status == SENT;}
  bool awaiting(void) const {return status == AWAITING;}
  bool ready(void) const {return status == READY;}
  bool received(void) const  {return status == RECEIVED;}
  bool empty(void) const {return status == EMPTY;}
  uint32_t get_grip_ID(void) const {return grip_ID;}

  void set_grip_ID(uint32_t _grip_ID) {grip_ID = _grip_ID;}

  private:
  std::vector<T> object;
  uint32_t status;
  uint32_t grip_ID; //! Only used in mesh passing
  enum {EMPTY, READY, SENT, AWAITING, RECEIVED};

};

#endif // buffer_h_
//---------------------------------------------------------------------------//
// end of buffer.h
//---------------------------------------------------------------------------//
