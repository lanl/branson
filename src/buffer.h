//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   buffer.h
 * \author Alex Long
 * \date   December 15 2015
 * \brief  Class for simplifying data management in MPI messaging
 * \note   ***COPYRIGHT_GOES_HERE****
 */
//---------------------------------------------------------------------------//

#ifndef buffer_h_
#define buffer_h_

#include <vector>
#include <mpi.h>

template<class T>
class Buffer {

  public:
  //! Default constructor
  Buffer()
  : status(EMPTY),
    buffer_recv_size(0),
    rank(MPI_PROC_NULL)
  {}

  //! Destructor
  ~Buffer() {}

  //! Fill the underlying buffer data with input vector
  void fill (std::vector<T> _object) {
    object = _object;
    status = READY;
  }

  //! Return pointer to the underlying buffer data (needed in MPI call)
  void *get_buffer (void) {
    if (object.size()>0) return &object[0];
    else return &object;
  }

  //! Return reference to buffer's vector object
  std::vector<T>& get_object(void) {return object;}

  //! Resize internal buffer vector
  void resize(uint32_t new_size) {object.resize(new_size);}

  //! Clear contents of internal buffer vector
  void clear(void) {object.clear(); status=EMPTY;}

  //! Set status to empty (ready to fill)
  void reset(void) {status=EMPTY;}

  //! Set status to sent (can be tested)
  void set_sent(void) {status = SENT;}

  //! Set status to received (data can be removed from buffer)
  void set_received(void) {status = RECEIVED;}

  //! Set status to awaiting (can be tested)
  void set_awaiting(void) {status = AWAITING;}

  //! Check to see if the status is "SENT"
  bool sent(void) const {return status == SENT;}

  //! Check to see if the status is "AWAITING"
  bool awaiting(void) const {return status == AWAITING;}

  //! Check to see if the status is "READY"
  bool ready(void) const {return status == READY;}

  //! Check to see if the status is "RECEIVED"
  bool received(void) const  {return status == RECEIVED;}

  //! Check to see if buffer is "EMPTY" 
  bool empty(void) const {return status == EMPTY;}

  //! Return the grip IDs that were received by this buffer
  const std::vector<uint32_t>& get_grip_IDs(void) const {return grip_IDs;}

  //! Return the single grip ID of the buffer in CELL_PASS_RMA mode
  uint32_t get_grip_ID(void) const {return grip_IDs[0];}

  //! Return the rank which the buffer is sending to or receiving from
  int32_t get_rank(void) const {return rank;}

  //! Return the actual size of the received message
  uint32_t get_receive_size(void) const {return buffer_recv_size;}

  //! Set the single grip ID associated with the buffer (CELL_PASS_RMA mode)
  void set_grip_ID(const uint32_t _grip_ID) {
    grip_IDs = std::vector<uint32_t>(1,_grip_ID);
  }

  //! Set the multiple grip IDs assocaited with this buffer (CELL_PASS mode)
  void set_grip_IDs(std::vector<uint32_t> _grip_IDs) {grip_IDs = _grip_IDs;}

  //! Set the rank which the buffer is sending to or receiving from
  void set_rank(uint32_t _rank) {rank=_rank;}

  //! Set actual size of receive buffer
  void set_receive_size(uint32_t _recv_size) {buffer_recv_size = _recv_size;}

  private:
  uint32_t status; //!< Current status of the buffer

  //! Actual elements sent over MPI (buffer is generally oversized)
  uint32_t buffer_recv_size;

  //! Rank (source for receive, destination for send) used in mesh passing
  int32_t rank;

  //! Grips received or sent, used for convenience in mesh passing
  std::vector<uint32_t> grip_IDs;

  std::vector<T> object; //!< Where sent/received data is stored

  //! Buffer statuses, used in completion routine
  enum {EMPTY, READY, SENT, AWAITING, RECEIVED};
};

#endif // buffer_h_
//---------------------------------------------------------------------------//
// end of buffer.h
//---------------------------------------------------------------------------//
