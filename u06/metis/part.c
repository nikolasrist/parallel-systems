/*==============================================================================
  
   Purpose:    small example on how to use Metis 5.x
   Author:     Rudolf Berrendorf
               Computer Science Department
               Bonn-Rhein-Sieg University of Applied Sciences
	       53754 Sankt Augustin, Germany
               rudolf.berrendorf@h-brs.de
  
==============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <libFHBRS.h>

#include <metis.h>		       // Metis 5 interface!

// idx_t is the int-type in Metis

/*----------------------------------------------------------------------------*/
// just for debugging: print the internal of a a graph

static void
print_graph (idx_t n_vertex,	       // number of vertices
	     idx_t * vwgt,	       // vertex weights
	     idx_t * xadj,	       // graph structure I
	     idx_t * adjncy,	       // graph structure II
	     idx_t * adjwgt	       // edge weights
	)
{
  printf ("graph (%ld vertices, %ld edges):\n", (long) n_vertex, (long) xadj[n_vertex] / 2);

  // print vertices with adjacency list
  printf ("%6s (%10s) : %s\n", "vertex", "weight", "neighbors(edge weight)");
  for (idx_t vertex = 0; vertex < n_vertex; vertex++)
    {
      // print vertex number and vertex weight
      printf ("%6ld (%10ld) : ", (long) vertex, (long) vwgt[vertex]);

      // adjacent vertices. print vertex number and edge weight
      idx_t start_edge = xadj[vertex];
      idx_t end_edge = xadj[vertex + 1] - 1;
      for (idx_t neighbor = start_edge; neighbor <= end_edge; neighbor++)
	printf ("%4ld(%ld) ", (long) adjncy[neighbor], (long) adjwgt[neighbor]);
      printf ("\n");
    }
}


/*----------------------------------------------------------------------------*/
// just for debugging: print out partitioning information

static void
print_partition (idx_t n_vertex,       // number of vertices
		 idx_t * vwgt,	       // vertex weights
		 idx_t * xadj,	       // graph structure I
		 idx_t * adjncy,       // graph structure II
		 idx_t * adjwgt,       // edge weights
		 idx_t part[n_vertex], // partitioning vector
		 idx_t nparts,	       // number of partitions
		 idx_t edgecut	       // edge cut
	)
{
  // partition
  printf ("partition result:\n");
  
  // print vertex numbers
  printf ("   vertex   : ");
  for (idx_t i = 0; i < n_vertex; i++)
    printf ("%2ld ", (long) i);
  printf ("\n");

  // print mapping to a parition
  printf ("on partition: ");
  for (idx_t i = 0; i < n_vertex; i++)
    printf ("%2ld ", (long) part[i]);
  printf ("\n");

  // edge cut
  printf ("edge cut=%ld\n", (long) edgecut);

  // compute vertex weights + crossing edge weights for all partitions
  printf ("costs for partitions: \n");
  float costs[nparts];
  for (idx_t i = 0; i < nparts; i++)
    costs[i] = 0;

  for (idx_t i = 0; i < n_vertex; i++)
    {
      idx_t my_part = part[i];

      // add vertex cost to partition
      costs[my_part] += vwgt[i];

      // determine costs for edges crossing partition from our vertex
      idx_t start_edge, end_edge;
      start_edge = xadj[i];
      end_edge = xadj[i + 1] - 1;
      for (idx_t j = start_edge; j <= end_edge; j++)
	{
	  if (my_part != part[adjncy[j]])
	    costs[part[i]] += adjwgt[j];
	}
    }

  // determine parition with maximum costs
  float maxcost = 0.0;
  for (idx_t i = 0; i < nparts; i++)
    {
      printf ("%f ", costs[i]);
      if (costs[i] > maxcost)
	maxcost = costs[i];
    }

  // maximum for all partitions is total weight
  printf ("\ntotal cost = %f\n", maxcost);
}


/*----------------------------------------------------------------------------*/
// small test program

int
main (int argc, char **argv)
{
  // see Metis documentation (p.20 ff for Metis 5.x)

  /*----------------------------------------------------------------------------*/
  // example graph

  // 0 -- 3
  // |    |
  // 1    4
  // |  / |
  // 2    5

  // number of vertices
  idx_t n_vertex = 6;

  // size of xadj is n_vertex+1.
  // This array stores for every vertex i in xadj[i] the index
  // where the adjacent vertices are found in adjncy[].
  // Looking at the example graph:
  // Neighbors of vertex 0 are found starting at adjncy[0] and ending in adjncy[2-1],
  // i.e. one less than the starting point of the next vertex.
  // Neighbors of vertex 1 are found starting at adjncy[2] and ending in adjncy[4-1],
  // i.e. one less than the starting point of the next vertex.
  // etc.
  idx_t xadj[] = { 0, 2, 4, 6, 8, 11, 12 };

  // number of edges (see below: adjcny stores every edge twice as (u,v) and (v,u))
  idx_t n_edge = 6;

  // size of adjncy is 2 * n_edge (every edge appears twice in undirected graph)
  idx_t adjncy[] = {
    1, 3,	                       // neighbors of vertex 0
    0, 2,			       // neighbors of vertex 1
    1, 4,			       // neighbors of vertex 2
    0, 4,			       // neighbors of vertex 3
    2, 3, 5,			       // neighbors of vertex 4
    4				       // neighbors of vertex 5
  };

  // vertex weights (size of vwgt is n_vertex)
  idx_t vwgt[] = { 20, 30, 80, 60, 100, 90 };

  // edge weights (size of adjwgt is n_edge*2; every edge exists twice for both "directions")
  idx_t adjwgt[] = { 10, 20, 10, 10, 10, 30, 20, 10, 30, 10, 10, 10 };

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
  idx_t nparts = 2;
  // !!!!!!!!!! you have to change this !!!!!!!!!!!!



  /*--------------------------------------------------------------------------*/

  // print header and input graph
  printf ("partitioning graph with %ld vertices and %ld edges for %ld partitions\n", (long)n_vertex, (long)n_edge, (long)nparts);
  print_graph (n_vertex, vwgt, xadj, adjncy, adjwgt);


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
  printf ("partitioning time: %.9f s\n", t0);
  // print parition
  print_partition (n_vertex, vwgt, xadj, adjncy, adjwgt, part, nparts, edgecut);

  return 0;
}

/*============================================================================*
 *                             that's all folks                               *
 *============================================================================*/
