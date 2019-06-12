//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   message_counter.h
 * \author Alex Long
 * \date   July 20 2016
 * \brief  Struct that holds all network message counters
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//----------------------------------------------------------------------------//

#ifndef message_counter_h_
#define message_counter_h_

#include <vector>

struct Message_Counter {

public:
  //! Constructor
  Message_Counter()
      : n_sends_posted(0), n_sends_completed(0),
        n_receives_posted(0), n_receives_completed(0) {}

  //! Destructor
  ~Message_Counter() {}

  //! Reset the counter data to zero
  void reset_counters(void) {
    n_sends_posted = 0;
    n_sends_completed = 0;
    n_receives_posted = 0;
    n_receives_completed = 0;
  }

  uint32_t n_sends_posted;       //!< Number of sent messages posted
  uint32_t n_sends_completed;    //!< Number of sent messages completed
  uint32_t n_receives_posted;    //!< Number of received messages completed
  uint32_t n_receives_completed; //!< Number of received messages completed
};

#endif // message_counter_h_
//----------------------------------------------------------------------------//
// end of message_counter.h
//----------------------------------------------------------------------------//
