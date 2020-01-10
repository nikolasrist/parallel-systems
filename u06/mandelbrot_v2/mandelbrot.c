/*==============================================================================
  
   Purpose:    MPI Mandelbrot (master/slave approach)
   Author:     Rudolf Berrendorf
               Computer Science Department
               Bonn-Rhein-Sieg University of Applied Sciences
	       53754 Sankt Augustin, Germany
               rudolf.berrendorf@h-brs.de
  
==============================================================================*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>

#include <mpi.h>

#include <libFHBRS.h>


//==============================================================================
// constants

#define NROW 1500                      // number of rows

#define	X_RESOLUTION NROW	       // number of pixels in x-direction (rows)
#define	Y_RESOLUTION NROW	       // number of pixels in y-direction (columns)

// debug output (0=nothing, 1=moderate, 2=more details)
#define DEBUG 0


//==============================================================================
// global variables

// MPI size and rank
static int size, rank;

// checksum if we don't display
static unsigned long checksum;
static double taskTimes[NROW];

//==============================================================================
/*  Only executed on master processor (�processor 0).
    Receive one row from a slave processor.
*/

static void receive_data () {
  int anziter[Y_RESOLUTION + 1];
  MPI_Status status;
  int err;


  // processor 0 receives information
  err = MPI_Recv (anziter, Y_RESOLUTION+1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
  assert(err == MPI_SUCCESS);

  // compute checksum (without line number in first place)
  for (int j = 1; j <= Y_RESOLUTION; j++) {
    checksum += anziter[j];
  }
}

static void receive_taskTime() {
  double taskTime[2];
  MPI_Status status;
  int err;
  // processor 0 receives task times
  err = MPI_Recv (taskTime, 2, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
  assert(err == MPI_SUCCESS);

  // printf("Received TaskTime: %.2f,%.2f\n", taskTime[0], taskTime[1]);
  // write task time in overall tasktime array
  taskTimes[(int)taskTime[0]] = taskTime[1];
}

//==============================================================================
/* Only executed on processors other than 0.
   Get gather data on a point (and in a non-batch environment we would display this pixel).
*/

static void drawPoint (int i, int j, int anziter) {
  // static: lifetime over all calls to the function
  static int iter[Y_RESOLUTION+1];

  /* processor 1..n gather information for a complete row
     and send the complete row to processor 0
  */

  // store row number
  iter[0] = i;

  // store value of (i,j)
  iter[j+1] = anziter;

  // last pixel in row?
  if (j == Y_RESOLUTION-1) {
    // end of row reached: send whole row to processor 0
    // we assume here that a processor computes all pixels of a row in sequence
    int err = MPI_Send (iter, Y_RESOLUTION+1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);
  }
}

//==============================================================================
/* does a block distribution from start-end to p processors returning in
   local_start / local_end block boundaries for processor iam
 */

static void
block_distribution (int start,	       /* start iteration */
		    int end,	       /* end iteration (incl.) */
		    int p,	       /* number of processors */
		    int iam,	       /* my processor number */
		    int *local_start,  /* local start iteration */
		    int *local_end     /* local end iteration (incl.) */
	)
{
  int n = end - start + 1;
  int q = n / p;
  int r = n - (q * p);

  if (r == 0) {
    *local_start = iam * q;
    *local_end = (iam + 1) * q - 1;
  } else {
    if (iam < r) {
      *local_start = iam * (q + 1);
      *local_end = *local_start + (q + 1) - 1;
    } else {
      *local_start = (r * (q + 1)) + (iam - r) * q;
      *local_end = *local_start + q - 1;
    }
  }

#if (DEBUG > 0)
  printf ("p=%d, n=%d, q=%d, r=%d, start=%d, end=%d, localstart=%d, localend=%d\n",
	  p, n, q, r, start, end, *local_start, *local_end);
#endif
}

//==============================================================================
// mandelbrot computation (on clients)

static void mandelbrot_client(int maxiter, double dx, double dy, double xmin, double ymin) {
  int start_iter, end_iter;

  /* all other processors compute rows.
     Determine start and end iteration for this processor
     by a simple block distribution in the reference implementation.

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     It makes sense to distribute the work to a process at this single
     place. Therefore replace the call to block_distribution with:
     1) get information which iterations this processors should work on
     i.e. get schedule vector from graph partitioning
     2) let this processor do all these assigned iteration of the outer
     loop
     i.e. make a new surrounding loop and check at the i-th place of the
     schedule vector whether this processor should execute iteration i (in
     this case: start_iter=end_iter=i)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  */
  block_distribution (0, X_RESOLUTION - 1, size - 1, rank - 1, &start_iter, &end_iter);
  
#if (DEBUG > 0)
  printf("process %d works from %d to %d\n", rank, start_iter, end_iter);
#endif
  

  // calculate values for every point in complex plane
  for (int i = start_iter; i <= end_iter; i++) {
    // measure row computation time, i.e. execution time for one single task
    double t_task = gettime();

#if (DEBUG > 1)
    printf("process %d starts working on row %d\n", rank, i);
#endif
  
    for (int j = 0; j < Y_RESOLUTION; j++) {
      int k;
      double absvalue, temp;
      struct {
        double real, imag;
      } z, c;

      // map point to window
      c.real = xmin + i * dx;
      c.imag = ymin + j * dy;
      z.real = z.imag = 0.0;
      k = 0;
      
      do {
        temp = z.real * z.real - z.imag * z.imag + c.real;
        z.imag = 2.0 * z.real * z.imag + c.imag;
        z.real = temp;
        absvalue = z.real * z.real + z.imag * z.imag;
        k++;
      } while (absvalue < 4.0 && k < maxiter);
      
      // display result (in our case just add to checksum)
      drawPoint (i, j, k);
    }

    // task time
    t_task = gettime() - t_task;
    double currentTask[2];
    currentTask[0] = i;
    currentTask[1] = t_task;
    int err = MPI_Send (currentTask, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);
  }
}


//==============================================================================
/* mandelbrot computation */

static void mandelbrot(int maxiter, double dx, double dy, double xmin, double ymin) {
  if (rank == 0) {
    // master processor waits for data from slave processors
    for (int i = 0; i < X_RESOLUTION; i++){
      receive_data ();
      receive_taskTime();
    }
  } else {
    // all clients work on mandelbrot computations and send results to master
    mandelbrot_client(maxiter, dx, dy, xmin, ymin);
  }
}


//==============================================================================
// main program

int
main (int argc, char **argv) {
  int maxiter, err;
  double xmin, ymin, xmax, ymax;
  double dx, dy;
  double t_start, t_end;


  // initialize MPI
  err = MPI_Init (&argc, &argv);
  assert(err == MPI_SUCCESS);
  err = MPI_Comm_size (MPI_COMM_WORLD, &size);
  assert(err == MPI_SUCCESS);
  err = MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  assert(err == MPI_SUCCESS);

  // initialization of mandelbrot variables
  xmin = -1.5;
  ymin = -1.5;
  xmax = 1.5;
  ymax = 1.5;
  maxiter = 100000;

  // init display variables
  dx = (xmax - xmin) / X_RESOLUTION;
  dy = (ymax - ymin) / Y_RESOLUTION;

  //--------------------------------------------------------------------------
  // mandelbrot computation

  // get start time
  t_start = gettime ();

  mandelbrot(maxiter, dx, dy, xmin, ymin);

  // get end time
  t_end = gettime ();

  //--------------------------------------------------------------------------

  if (rank == 0) {
    printf ("calculation took %.2f s on %d+1 processors\n", t_end - t_start, size-1);
    printf("checksum = %lu\n", checksum);
    printf("task times: ");
    for(int i=0; i < NROW; i++) {
      printf("row %d took %.2f\n", i, taskTimes[i]);
    }
  }

  // exit MPI
  err = MPI_Finalize ();
  assert(err == MPI_SUCCESS);

  return EXIT_SUCCESS;
}

/*============================================================================*
 *                             that's all folks                               *
 *============================================================================*/
