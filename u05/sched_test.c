/*==============================================================================
  
   Purpose          : test program for loop scheduling algorithms
   Author           : Rudolf Berrendorf
                      Computer Science Department
                      Bonn-Rhein-Sieg University of Applied Sciences
	              53754 Sankt Augustin, Germany
                      rudolf.berrendorf@h-brs.de
  
==============================================================================*/
	 
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <libFHBRS.h>

#if defined(_OPENMP)
#include <omp.h>
#endif
	
#include "sched.h"
	
long int random(void);
void srandom(unsigned int seed);

/*============================================================================*/ 
/* macros */ 
	
#define N_WORK 4                        /* number of work functions */
#define N_SCHED 6                       /* number of scheduling algorithms */

/* computational work inside work functions */
#define calc(i) ((i % 10) + (i % 2))


/*============================================================================*/ 
/* types */ 
	
/* C type of a work function */ 
typedef void (*work_t) (int, int *);


/*============================================================================*/ 
/* variables */ 

/* a profile consists of a pair number of iterations to execute and work load
   of an iteration */
/* number of iterations to be executed (-1 = end of list) */
static int n_iterations[] = { 1000, 100000, 10000000 , -1};

/* workload factor (smaller = less workload) */
static float work_factors[] = { 1000.0, 0.1, 0.00001 , 0.0 };

/* global variable with number of loop iterations */
static int n;

/* global variable that gets one of the work_factor values assigned */
static float work_factor;

/* we compute some value in this variable in all working functions */
static volatile int global_x = 0;

/* some random numbers between 0.5 and 1.1 for random working strategy */
static float *random_value;


/* scheduling setup functions (called once before every loop scheduling) */
static scheduling_setup_function_t *setup[N_SCHED] = 
        { sequential_scheduling_setup,
	  block_scheduling_setup,
	  blockcyclic_scheduling_setup,
	  self_scheduling_setup,
	  gss_setup,
	  factoring_setup
};

/* scheduling functions */
static scheduling_function_t *sched[N_SCHED] = 
	{ sequential_scheduling,
	  block_scheduling,
	  blockcyclic_scheduling,
	  self_scheduling,
	  gss,
	  factoring
};

/* names of scheduling function */
static char *sched_names[N_SCHED] = 
	{ "sequential",
	  "block",
	  "blockcyclic",
	  "self",
	  "gss",
	  "factoring"
};


/*============================================================================*/ 
/* work load function: independent of iteration number produce same work load */

static void
work_constant (int iter, int *result) 
{
  for (int i = 0; i < (int)((n/2) * work_factor); i++)
    *result += calc(i);
} 

/*============================================================================*/ 
/*  work load function: dependent of iteration number produce increasing work load */

static void
work_increasing (int iter, int *result) 
{
  for (int i = 0; i < (int)(iter * work_factor); i++)
    *result += calc(i);
} 

/*============================================================================*/ 
/*  work load function: dependent of iteration number produce decreasing work load */

static void
work_decreasing (int iter, int *result) 
{
  for (int i = 0; i < (int)((n - iter + 1) * work_factor); i++)
    *result += calc(i);
} 

/*============================================================================*/ 
/*  work load function: dependent of iteration number produce random work load */

static void
work_random (int iter, int *result) 
{
  for (int i = 0; i < (int)((n/2) * (random_value[iter] * work_factor)); i++)
    *result += calc(i);
} 

/*============================================================================*/ 

void dummy(int *g)
{
}

/*============================================================================*/ 

int
main (int argc, char **argv) 
{
  /* reference result values to check for correctness */
  int correct_result[N_WORK];
  /* work functions */
  work_t work_funs[N_WORK] = {work_constant, work_increasing, work_decreasing, work_random};


  /* dummy call for initialization */
  (void) gettime ();

  
  /*------------------------------------------------------------------------*/
  /* test loop: we iterate over:
     - work profiles (this loop)
     - work functions (next inner loop)
     - scheduling algorithms (innermost loop)
  */

  for(int niter = 0; n_iterations[niter] != -1; niter++)
    {

      /* test settings */
      /* number of iterations */
      n = n_iterations[niter];
      /* work load factor for work done in iterations */
      work_factor = work_factors[niter];

      printf("%d iterations, %f workload factor\n", n, work_factor);


      /*------------------------------------------------------------------------*/
      /* produce some random numbers between 0.5 and 1.5 for random working strategy */

      random_value = malloc(n * sizeof(*random_value));
      srandom(0);
      for(int i=0; i<n; i++)
	random_value[i] = (random() / (float)RAND_MAX) + 0.5f;


      /*------------------------------------------------------------------------*/
      /* produce correct reference results for all working functions */ 

      for(int work_fun=0; work_fun < N_WORK; work_fun++)
	{
	  global_x = 0;
	  for (int i = 0; i < n; i++)
	    (work_funs[work_fun]) (i, (int *)&global_x);
	  /* store rerefence value */
	  correct_result[work_fun] = global_x;
	}
      
      
      /*------------------------------------------------------------------------*/
      /* here starts the parallel region */

#pragma omp parallel
      {
#if defined(_OPENMP)
	int iam = omp_get_thread_num();
	int p = omp_get_num_threads();
#else
	/* one processor */
	int iam = 0;
	int p = 1;
#endif

	/* thread private variable for partial result */
	int my_global_x;


#pragma omp master
	  printf("%d processor(s), times in ms\n%-15s %10s %10s %10s %10s\n",
		 p, "sched-alg", "constant", "increasing", "decreasing","random");

	/*----------------------------------------------------------------------*/
	/* iterate over scheduling algorithms */

	for (int alg = 0; alg < N_SCHED; alg++)  
	  {
#pragma omp master
	    printf("%-15s", sched_names[alg]);
	    
	    /*------------------------------------------------------------------*/
	    /* iterate over work loads */

	    for(int work_fun=0; work_fun < N_WORK; work_fun++)
	      {
		int start_iteration = INVALID_ITERATION, end_iteration = INVALID_ITERATION;
		
#pragma omp master
		{
		  /* initialize dummy variable to be computed by work function */ 
		  global_x = 0;
		  
		  /* specific scheduling algorithm setup */ 
		  setup[alg](n, p, iam);
		}
		
#pragma omp barrier
		/* on every processor */
		my_global_x = 0;
		
		/* start time */ 
		double t0 = gettime ();
		
		/* run */ 
		while (sched[alg] (&start_iteration, &end_iteration, n, p, iam))
		  {
#if defined(DEBUG)
		    printf("%d-%d ", start_iteration, end_iteration);
		    fflush(stdout);
#endif
		    
		    /* do work on chunk assigned by scheduliong algorithm */ 
		    for (int i = start_iteration; i <= end_iteration; i++)
		      /* do the work */ 
		      work_funs[work_fun] (i, &my_global_x);
		  } 
		
#pragma omp barrier
		/* wait untiall processors have finished */
#pragma omp master
		{
		  /* end time */ 
		  t0 = gettime () - t0;
		  printf(" %10.3f", t0 * 1e3);
		}
		
#pragma omp critical
		global_x += my_global_x;
#pragma omp barrier
		;
#pragma omp master
		{
		  /* check for correctness */ 
		  if (global_x != correct_result[work_fun])
		    printf ("error: %d != %d\n", global_x, correct_result[work_fun]);
		}
		dummy((int *)&global_x);
		
	      } /* for(int work_fun=0;...) */
	    
#pragma omp master
	    printf("\n");
	    
	  } /* for(int alg...) */
      } /* pragma omp parallel */

      free(random_value);

    } /* for(int niter...) */
      

  return 0;
}
  
  
/*============================================================================*
 *                             that's all folks                               *
 *============================================================================*/ 
