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

#include <metis.h>		       // Metis 5 interface!

// constants

#define NROW 1500                     // number of rows

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
static unsigned long checksum_v2;
static double taskTimes[X_RESOLUTION];
static int partitions[X_RESOLUTION];

//==============================================================================
// metis helper
//==============================================================================


// static void
// print_graph (idx_t n_vertex,	       // number of vertices
// 	     idx_t * vwgt,	       // vertex weights
// 	     idx_t * xadj,	       // graph structure I
// 	     idx_t * adjncy,	       // graph structure II
// 	     idx_t * adjwgt	       // edge weights
// 	)
// {
//   printf ("graph (%ld vertices, %ld edges):\n", (long) n_vertex, (long) xadj[n_vertex] / 2);

//   // print vertices with adjacency list
//   printf ("%6s (%10s) : %s\n", "vertex", "weight", "neighbors(edge weight)");
//   for (idx_t vertex = 0; vertex < n_vertex; vertex++)
//     {
//       // print vertex number and vertex weight
//       printf ("%6ld (%10ld) : ", (long) vertex, (long) vwgt[vertex]);

//       // adjacent vertices. print vertex number and edge weight
//       idx_t start_edge = xadj[vertex];
//       idx_t end_edge = xadj[vertex + 1] - 1;
//       for (idx_t neighbor = start_edge; neighbor <= end_edge; neighbor++)
// 	printf ("%4ld(%ld) ", (long) adjncy[neighbor], (long) adjwgt[neighbor]);
//       printf ("\n");
//     }
// }


/*----------------------------------------------------------------------------*/
// just for debugging: print out partitioning information

// static void
// print_partition (idx_t n_vertex,       // number of vertices
// 		 idx_t * vwgt,	       // vertex weights
// 		 idx_t * xadj,	       // graph structure I
// 		 idx_t * adjncy,       // graph structure II
// 		 idx_t * adjwgt,       // edge weights
// 		 idx_t part[n_vertex], // partitioning vector
// 		 idx_t nparts,	       // number of partitions
// 		 idx_t edgecut	       // edge cut
// 	)
// {
//   // partition
//   printf ("partition result:\n");
  
//   // print vertex numbers
//   printf ("   vertex   : ");
//   for (idx_t i = 0; i < n_vertex; i++)
//     printf ("%2ld ", (long) i);
//   printf ("\n");

//   // print mapping to a parition
//   printf ("on partition: ");
//   for (idx_t i = 0; i < n_vertex; i++)
//     printf ("%2ld ", (long) part[i]);
//   printf ("\n");

//   // edge cut
//   printf ("edge cut=%ld\n", (long) edgecut);

//   // compute vertex weights + crossing edge weights for all partitions
//   printf ("costs for partitions: \n");
//   float costs[nparts];
//   for (idx_t i = 0; i < nparts; i++)
//     costs[i] = 0;

//   for (idx_t i = 0; i < n_vertex; i++)
//     {
//       idx_t my_part = part[i];

//       // add vertex cost to partition
//       costs[my_part] += vwgt[i];

//       // determine costs for edges crossing partition from our vertex
//       idx_t start_edge, end_edge;
//       start_edge = xadj[i];
//       end_edge = xadj[i + 1] - 1;
//       for (idx_t j = start_edge; j <= end_edge; j++)
// 	{
// 	  if (my_part != part[adjncy[j]])
// 	    costs[part[i]] += adjwgt[j];
// 	}
//     }

//   // determine parition with maximum costs
//   float maxcost = 0.0;
//   for (idx_t i = 0; i < nparts; i++)
//     {
//       printf ("%f ", costs[i]);
//       if (costs[i] > maxcost)
// 	    maxcost = costs[i];
//     }

//   // maximum for all partitions is total weight
//   printf ("\ntotal cost = %f\n", maxcost);
// }

static
void calculatePartitions(double taskTimes[]){

  // number of vertices
  idx_t n_vertex = X_RESOLUTION;

  // size of xadj is n_vertex+1.
  // This array stores for every vertex i in xadj[i] the index
  // where the adjacent vertices are found in adjncy[].
  // Looking at the example graph:
  // Neighbors of vertex 0 are found starting at adjncy[0] and ending in adjncy[2-1],
  // i.e. one less than the starting point of the next vertex.
  // Neighbors of vertex 1 are found starting at adjncy[2] and ending in adjncy[4-1],
  // i.e. one less than the starting point of the next vertex.
  // etc.

  // calculate one adjacency for each node
  idx_t xadj[n_vertex+1];
  idx_t edgeIndex = 0;
  for(idx_t i = 0; i<=n_vertex; i++){
    if(i == 0) {
      xadj[i] = edgeIndex;
      edgeIndex++;  
    }else if (i == n_vertex){
      xadj[i] = edgeIndex-1;
    } else {
      xadj[i] = edgeIndex;
      edgeIndex += 2;
    }
  }
  
  idx_t n_edge = X_RESOLUTION - 1;

  // size of adjncy is 2 * n_edge (every edge appears twice in undirected graph)
  idx_t adjncy_size = (n_edge * 2);
  idx_t adjncy[adjncy_size];

  int writeIndex = 0;
  for(idx_t i = 0; i < n_vertex; i++) {
    if(i == 0) { // first node only has i+1 as neighbour
      // printf("First node: %d\n", 1);
      adjncy[0] = 1;
      writeIndex++;
    } else if (i == n_vertex - 1) {
      // last node only has i-1 as neighbour
      // printf("Last node: %d\n", n_vertex-1);
      adjncy[adjncy_size-1] = n_vertex-2;
    } else {
      // inner nodes have two neighbours i+1 and i-1
      // printf("Inner Node - calculated node %d - written: %d, %d to index %d\n", i, i-1, i+1, writeIndex);
      adjncy[writeIndex] = i-1;
      adjncy[writeIndex+1] = i+1;
      writeIndex += 2;
    }
    
  }

  idx_t vwgt[n_vertex];

  for(idx_t i = 0; i < n_vertex; i++ ){
    idx_t tmp = (idx_t)(taskTimes[i] * 1000);
    if(tmp == 0){
      tmp = 1;
    }
    vwgt[i] = tmp;
  }

  // edge weights (size of adjwgt is n_edge*2; every edge exists twice for both "directions")
  idx_t adjwgt[n_edge*2];

  for(idx_t i=0; i<n_edge*2;i++){
    adjwgt[i] = 1;
  }

  // the resulting partition:
  // for every vertex i we find in part[i] afterwards the partition number for that vertex
  idx_t part[n_vertex];


  // Metis options (initialized to default values)
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions (options);
  // partitioning method: k-way partitioning
  options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
  // edge cut minimization
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  // C-style numbering
  options[METIS_OPTION_NUMBERING] = 0;

  // number of balancing constraints
  idx_t ncon = 1;

  // edge cut after partitioning
  idx_t edgecut;
  // !!!!!!!!!! you have to change this !!!!!!!!!!!!
  // number of partitions we want to partition to
  idx_t nparts = size - 1;
  // !!!!!!!!!! you have to change this !!!!!!!!!!!!
  /*--------------------------------------------------------------------------*/
  // print header and input graph
  // printf ("partitioning graph with %ld vertices and %ld edges for %ld partitions\n", (long)n_vertex, (long)n_edge, (long)nparts);
  // print_graph (n_vertex, vwgt, xadj, adjncy, adjwgt);

  // call graph partitioning
  double t0 = gettime ();
  int rc = METIS_PartGraphKway (&n_vertex,	// number of vertices
				&ncon,          // number of balancing constraints
				xadj,           // adjacency structure of graph
				adjncy,	        // adjacency structure of graph
				vwgt,           // vertex weights
				NULL,           // only for total communication volume
				adjwgt,	        // edge weights
				&nparts,	// number of partitions wanted
				NULL,           // tpwgts: no desired partition weights
				NULL,           // ubvec: allowed load imbalance
				options,	// special options
				&edgecut,	// objective value (edge cut)
				part     // vector with partition information for each vertex
	  );
  t0 = gettime () - t0;

  // check Metis return value
  switch (rc)
    {
      case METIS_OK:
	break;
      case METIS_ERROR_INPUT:
	printf ("error in Metis input\n");
	exit (1);
      case METIS_ERROR_MEMORY:
	printf ("no more memory in Metis\n");
	exit (1);
      case METIS_ERROR:
	printf ("some error in Metis\n");
	exit (1);
      default:
	printf ("unknown return code ffrom Metis\n");
	exit (1);
    }


  /*--------------------------------------------------------------------------*/
  /* print results */

  // time for partitioning
  // printf ("partitioning time: %.9f s\n", t0);
  // print parition
  // print_partition (n_vertex, vwgt, xadj, adjncy, adjwgt, part, nparts, edgecut);

  // printf("Write partitions into partitions array.\n");
  for(int i = 0; i < n_vertex; i++) {
    // printf("Write Partition: %d\n", part[i]);
    partitions[i] = part[i] + 1;
    // printf("Saved Partition: %d\n", partitions[i]);
  }
  // printf("Done.\n");
} 

//==============================================================================
//==============================================================================
// NORMAL PROGRAM START POINT
//==============================================================================
//==============================================================================


//==============================================================================
/*  Only executed on master processor (processor 0).
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

static void receive_data_v2 () {
  int anziter[Y_RESOLUTION + 1];
  MPI_Status status;
  int err;
  // processor 0 receives information
  err = MPI_Recv (anziter, Y_RESOLUTION+1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
  assert(err == MPI_SUCCESS);

  // compute checksum (without line number in first place)
  for (int j = 1; j <= Y_RESOLUTION; j++) {
    checksum_v2 += anziter[j];
  }
}

static void receive_taskTime() {
  double taskTime[2];
  MPI_Status status;
  int err;
  // processor 0 receives task times
  err = MPI_Recv (taskTime, 2, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
  // printf("ReceiveTaskTime - err=%d\n", err);
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

// static void
// block_distribution (int start,	       /* start iteration */
// 		    int end,	       /* end iteration (incl.) */
// 		    int p,	       /* number of processors */
// 		    int iam,	       /* my processor number */
// 		    int *local_start,  /* local start iteration */
// 		    int *local_end     /* local end iteration (incl.) */
// 	)
// {
//   int n = end - start + 1;
//   int q = n / p;
//   int r = n - (q * p);

//   if (r == 0) {
//     *local_start = iam * q;
//     *local_end = (iam + 1) * q - 1;
//   } else {
//     if (iam < r) {
//       *local_start = iam * (q + 1);
//       *local_end = *local_start + (q + 1) - 1;
//     } else {
//       *local_start = (r * (q + 1)) + (iam - r) * q;
//       *local_end = *local_start + q - 1;
//     }
//   }

// #if (DEBUG > 0)
//   printf ("p=%d, n=%d, q=%d, r=%d, start=%d, end=%d, localstart=%d, localend=%d\n",
// 	  p, n, q, r, start, end, *local_start, *local_end);
// #endif
// }

//==============================================================================
// mandelbrot computation (on clients)

static void mandelbrot_client(int maxiter, double dx, double dy, double xmin, double ymin) {
  // int start_iter, end_iter;

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
  // block_distribution (0, X_RESOLUTION - 1, size - 1, rank - 1, &start_iter, &end_iter);
  
#if (DEBUG > 0)
  printf("process %d works from %d to %d\n", rank, start_iter, end_iter);
#endif
  

  // calculate values for every point in complex plane
  for (int i = 0; i < X_RESOLUTION; i++) {
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

    // send task time to processor
    double currentTask[2];
    currentTask[0] = i;
    currentTask[1] = t_task;
    // printf("Processor %d sends taskTime: %.5f of task: %d\n",rank, t_task, i);
    int err = MPI_Send (currentTask, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);
  }
}

//==============================================================================
// mandelbrot v2 computation (on clients)

static void mandelbrot_client_v2(int maxiter, double dx, double dy, double xmin, double ymin) {
  // int start_iter, end_iter;

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
  // block_distribution (0, X_RESOLUTION - 1, size - 1, rank - 1, &start_iter, &end_iter);
  
#if (DEBUG > 0)
  printf("process %d works from %d to %d\n", rank, start_iter, end_iter);
#endif
  

  // calculate values for every point in complex plane
  for (int i = 0; i < X_RESOLUTION; i++) {
    if(partitions[i] != rank) {
      // printf("Row is not in my partition: part=%d, rank=%d\n", partitions[i], rank);
      continue;
    }
    // measure row computation time, i.e. execution time for one single task
    double t_task = gettime();
    // printf("process %d checks working on row %d with partition=%d\n", rank, i, partitions[i]);
    // printf("process %d starts working on row %d\n", rank, i);
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
  }
}

//==============================================================================
/* mandelbrot computation */

static void mandelbrot(int maxiter, double dx, double dy, double xmin, double ymin) {
  if (rank == 0) {
    // master processor waits for data from slave processors
    // printf("Client:%d receiving tasktimes\n", rank);
    for (int i = 0; i < X_RESOLUTION; i++){
      // printf("Client:%d receiving tasktimes - row:%d\n", rank, i);
      receive_data ();
      receive_taskTime();
    }
  } else if (rank == 1) {
    // client 1 work on mandelbrot computations and send results to master
    // printf("Client:%d, calulating task times.\n", rank);
    mandelbrot_client(maxiter, dx, dy, xmin, ymin);
  }
}

//==============================================================================
/* mandelbrot computation */

static void mandelbrot_v2(int maxiter, double dx, double dy, double xmin, double ymin) {
  // printf("Starting Mandelbrot_v2 - %d\n", rank);
  if (rank == 0) {
    // printf("Master is going to receive all the data! %d\n", rank);
    // master processor waits for data from slave processors
    for (int i = 0; i < X_RESOLUTION; i++){
      // printf("Receive data for row: %d\n", i);
      receive_data_v2 ();
    }
  } else {
    // all clients work on mandelbrot computations and send results to master
    // printf("Call client function - %d\n", rank);
    mandelbrot_client_v2(maxiter, dx, dy, xmin, ymin);
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
    // printf("Master triggers partitioning\n");
    if(size > 2){
      // printf("Calculate paritions with metis\n");
      calculatePartitions(taskTimes);
    } else {
      // printf("Only sequential working - skip metis\n");
      for(int i = 0; i < X_RESOLUTION; i++){
        partitions[i] = 1;
      }
    }
    for(int i = 1; i < size; i++) {
      // printf("Send partitions to %d", i);
      int err = MPI_Send (partitions, X_RESOLUTION, MPI_INT, i, 0, MPI_COMM_WORLD);
      assert(err == MPI_SUCCESS);
    }
    // printf("Master ist done.\n");
  }
  if(rank != 0){
    MPI_Status status;
    // printf("Receive partitions from master. %d", rank);
    err = MPI_Recv (partitions, X_RESOLUTION, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
    assert(err == MPI_SUCCESS);
  }
  //--------------------------------------------------------------------------
  // mandelbrot computation with graph partioning
  // printf("Proc: %d waiting for barrier.\n", rank);
  MPI_Barrier(MPI_COMM_WORLD);

  // printf("Start new run with partitioning. %d\n", rank);
  // get start time
  t_start = gettime ();

  mandelbrot_v2(maxiter, dx, dy, xmin, ymin);

  // get end time
  t_end = gettime ();

  //--------------------------------------------------------------------------
  // mandelbrot computation with graph partioning
  if (rank == 0) {
    printf ("calculation v2 took %.2f s on %d+1 processors\n", t_end - t_start, size-1);
    printf("checksum = %lu\n", checksum_v2);
    // printf("task times: ");
    // for(int i=0; i < NROW; i++) {
    //   printf("row %d took %.2f\n", i, taskTimes[i]);
    // }
  }

  //--------------------------------------------------------------------------

    // exit MPI
  err = MPI_Finalize ();
  assert(err == MPI_SUCCESS);

  return EXIT_SUCCESS;
}

/*============================================================================*
 *                             that's all folks                               *
 *============================================================================*/
