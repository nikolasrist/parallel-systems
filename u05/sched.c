/*==============================================================================
  
   Purpose          : loop scheduling algorithms
   Author           : Rudolf Berrendorf
                      Computer Science Department
                      Bonn-Rhein-Sieg University of Applied Sciences
	              53754 Sankt Augustin, Germany
                      rudolf.berrendorf@h-brs.de
  
==============================================================================*/
	
#include <pthread.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <limits.h>

#if defined(_OPENMP)
#include <omp.h>
#endif
	
#include "sched.h"
	

/*============================================================================*/ 
/*
  Some remarks:
  -  every scheduling algorithm has an initialization function (xx_setup)
     and a scheduling function (xxx).
     The initialization gets called once, the scheduling function
     possibly several times and for several parallel loops

  -  a scheduling function has the prototype:
     bool xxx (int *start_iteration, int *end_iteration, int n, int p, int iam) 
     where *start_iteration and *end_iteration  must be assigned meaningful
     values for the iteration block assigned to the calling processor
     n is the total number of iteration
     p the number of threads
     iam my thread number between 0 and p-1

  -  a scheduling function is called every first time for one loop scheduling
     with *start_iteration==INVALID_ITERATION. I.e. you can see the start of a
     new loop to be scheduled if *start_iteration==INVALID_ITERATION.
     You may use that e.g. to reset some variables etc.
     See for example sequential_scheduling.

  -  *start_iteration and *end_iteration have the same values when the
     scheduling function is called again (other than the first time as
     described above)

  -  A scheduling function must *assign start_iteration and *end_iteration
     values for the next iteration block to be executed.
     E.g. block scheduling of 100 iterations on 2 processors would return
     0-49 when called by the first processor.

  -  A scheduling function returns true if start_iteration/end_iteration have
     meaningful values, i.e. you have assigned another chunk to that processor.
     The function returns false if no iteration is left over for the calling
     processor.

  -  Loops to be scheduled run always from 0 to some n

  -  There is a static variable current_iteration that you may use for dynamic
     scheduling algorithms e.g. to keep track of the last scheduled iteration

  Example call flow:
    xxx_setup()

    // new parallel loop
    xxx(*start_iteration==INVALID_ITERATION, *end_iteration, ...)      -> true
    xxx(*start_iteration==some_value, *end_iteration==some_value, ...) -> true
    xxx(*start_iteration==some_value, *end_iteration==some_value, ...) -> false

    // new parallel loop
    xxx(*start_iteration==INVALID_ITERATION, *end_iteration, ...)      -> true
    xxx(*start_iteration==some_value, *end_iteration==some_value, ...) -> false
 */

/*============================================================================*/ 
/* helpful macros */ 
	
#define min(x,y) (((x)<(y))?(x):(y))
#define max(x,y) (((x)>(y))?(x):(y))	

/*============================================================================*/ 
/* global variables */ 


/* current iteration (you need that for dynamic algorithms
   to know what was scheduled last) */
static volatile int current_iteration;


/*============================================================================*/ 

static int upper_gauss(int x, int y) {
  int div = x / y;

  if(div * y == x)
    return div;
  else
    return div + 1;
}

/*============================================================================*/ 
/* sequential (all iterations get assigned to the first calling processor) */ 

void
sequential_scheduling_setup (int n, int p, int iam) {
  /* nothing */
}

bool sequential_scheduling (int *start_iteration, int *end_iteration, int n, int p, int iam)  {
  bool first_time = (*start_iteration == INVALID_ITERATION);

  if (first_time && (iam == 0)) {
    /* we schedule only once all iterations and only on first processor */ 
    *start_iteration = 0;
    *end_iteration = n-1;
    return true;

  } else {
    return false;
  }
}


/*============================================================================*/ 
/* static block scheduling */ 

static volatile int block_q;
static volatile int block_r;

void
block_scheduling_setup (int n, int p, int iam) {
  /* block factor (iterations from 0 to n-1) */ 
  block_q = (n - 1) / p;
  
  /* number of processors that get blocks of size (q+1) */ 
  block_r = n - (block_q * p);
}

bool block_scheduling (int *start_iteration, int *end_iteration, int n, int p, int iam)  {
  bool first_time = (*start_iteration == INVALID_ITERATION);

  if (iam <= block_r-1) {
    /* first r processors get larger blocks of size q+1 */ 
    *start_iteration = iam * (block_q + 1);
    *end_iteration = ((iam + 1) * (block_q + 1)) - 1;
  } else {
    /* rest gets smaller blocks of size q */ 
    *start_iteration = block_r * (block_q + 1)
      + ((iam - block_r) * block_q);
    *end_iteration = block_r * (block_q + 1)
      + ((iam - block_r + 1) * block_q) - 1;
  }
  
  if (first_time) {
    /* we schedule only once */ 
    return true;
  } else {
    return false;
  }
}


/*============================================================================*/ 
/* static block-cyclic scheduling */ 

static volatile int cyclic_q;

void
blockcyclic_scheduling_setup (int n, int p, int iam) {
  /* block factor constantly 1 */ 
  cyclic_q = 1;
}

bool blockcyclic_scheduling (int *start_iteration, int *end_iteration, int n, int p, int iam) {  
  if(*start_iteration == INVALID_ITERATION) {
    *start_iteration = iam * cyclic_q;
    *end_iteration = min(*start_iteration + cyclic_q - 1, n-1);
  } else {
    *start_iteration += cyclic_q * p;
    *end_iteration = min(*end_iteration + cyclic_q * p, n-1);
  }

  return *start_iteration < n;
}


/*============================================================================*/ 
/* dynamic self scheduling */ 

void self_scheduling_setup (int n, int p, int iam) {
  current_iteration = 0;
  /* your code comes here */
}

bool self_scheduling (int *start_iteration, int *end_iteration, int n, int p, int iam) {
    bool nextIter = false;
    #pragma omp critical
    {
      if(current_iteration < n){
        *start_iteration = current_iteration;
        *end_iteration = current_iteration++;
        nextIter = true;
      } 
    }
  return nextIter;
}


/*============================================================================*/ 
/* dynamic guided self scheduling */ 

static volatile int remainingIters;
static volatile int chunkSize;

void
gss_setup (int n, int p, int iam) {
  current_iteration = 0;
  remainingIters = n;
  chunkSize = 0;

  /* your code comes here */
}

bool gss (int *start_iteration, int *end_iteration, int n, int p, int iam) {  
  /* your code comes here */
  bool hasNewIter = false;
  bool firstChunk = current_iteration == 0;
  #pragma omp critical
  {
    if(current_iteration < n) {
      chunkSize = upper_gauss(remainingIters, p);
      *start_iteration = current_iteration;
      current_iteration += chunkSize;
      *end_iteration = current_iteration - 1;
      remainingIters = remainingIters - chunkSize;
      // printf("ChunkSize: %d, RemainingIters: %d, StartPoint: %d, EndPoint: %d\n", chunkSize, remainingIters, *start_iteration, *end_iteration);
      hasNewIter = true;
    }
  }
  return hasNewIter;
}


/*============================================================================*/ 
/* dynamic factoring (simplified version 1) */ 

static volatile int fact_remaining_iters;
static volatile int fact_chunk_size;
static volatile int fact_chunkCalls;

void
factoring_setup (int n, int p, int iam) {
  current_iteration = 0;
  fact_chunkCalls = 0;
  fact_remaining_iters = n;
  fact_chunk_size = 0; 
  /* your code comes here */
}

bool factoring (int *start_iteration, int *end_iteration, int n, int p, int iam) {
  bool first_time = (*start_iteration == INVALID_ITERATION);
  bool hasNewIter = false;
  bool firstChunk = current_iteration == 0;
  #pragma omp critical
  { 
    if(fact_chunkCalls == 0 || fact_chunkCalls == p) {
      fact_chunk_size = upper_gauss(fact_remaining_iters, 2*p);
      fact_remaining_iters = fact_remaining_iters - (p*fact_chunk_size);
      fact_chunkCalls = 0;
    }

    if(current_iteration < n) {
      *start_iteration = current_iteration;
      current_iteration += fact_chunk_size;
      *end_iteration = current_iteration - 1;
      // printf("ChunkSize: %d, RemainingIters: %d, StartPoint: %d, EndPoint: %d\n", fact_chunk_size, fact_remaining_iters, *start_iteration, *end_iteration);
      fact_chunkCalls++;
      hasNewIter = true;
    }
  }
  return hasNewIter;
}


/*============================================================================*
 *                             that's all folks                               *
 *============================================================================*/ 
