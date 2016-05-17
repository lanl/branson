
#include <time.h>
#include <sys/time.h>

#ifndef timing_functions_h_
#define timing_functions_h_

void print_elapsed_inside(const char* desc, struct timeval* start, 
  struct timeval* end) 
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
  printf("%s  %f \n", desc, 
    (elapsed.tv_sec*1000000 + elapsed.tv_usec)/1000000.0 );
}


void print_elapsed_inside_w_rank(const char* desc, struct timeval* start, 
  struct timeval* end, const int& rank) 
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
  printf("rank :%d %s  %f \n", rank, desc, 
    (elapsed.tv_sec*1000000 + elapsed.tv_usec)/1000000.0 );
}

double get_runtime(struct timeval* start,  struct timeval* end) {
  struct timeval elapsed;

  // calculate elapsed time
  if(start->tv_usec > end->tv_usec) 
  {
    end->tv_usec += 1000000;
    end->tv_sec--;
  }
  elapsed.tv_usec = end->tv_usec - start->tv_usec;
  elapsed.tv_sec  = end->tv_sec  - start->tv_sec;
  return double((elapsed.tv_sec*1000000 + elapsed.tv_usec)/1000000.0);
}


#endif // timing_functions_h_
