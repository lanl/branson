
#include <time.h>
#include <sys/time.h>

#ifndef timing_functions_h_
#define timing_functions_h_

void print_elapsed_inside(const char* desc, struct timeval* start, struct timeval* end) 
{
  //prints timing statistics
  struct timeval elapsed;
  /* calculate elapsed time */
  if(start->tv_usec > end->tv_usec) 
  {
    end->tv_usec += 1000000;
    end->tv_sec--;
  }
  elapsed.tv_usec = end->tv_usec - start->tv_usec;
  elapsed.tv_sec  = end->tv_sec  - start->tv_sec;
  printf("%s  %f \n", desc, (elapsed.tv_sec*1000000 + elapsed.tv_usec)/1000000.0 );
}

#endif // timing_functions_h_
