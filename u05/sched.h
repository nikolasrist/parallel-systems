/*==============================================================================
  
   Purpose          : loop scheduling algorithms
   Author           : Rudolf Berrendorf
                      Computer Science Department
                      Bonn-Rhein-Sieg University of Applied Sciences
	              53754 Sankt Augustin, Germany
                      rudolf.berrendorf@h-brs.de
  
==============================================================================*/

#if !defined(SCHED_H)
#define SCHED_H

#include <limits.h>


/*============================================================================*/
/* macros */

/* some value indicating an invalid iteration number */
#define INVALID_ITERATION INT_MIN


/*============================================================================*/
/* types */

/* setup a scheduling setup function
   Such a function initializes everything you need for a scheduler function
 */
typedef void (scheduling_setup_function_t) (int n, int p, int iam);

/* schedule function
   A schedule function returns a chunk of work in form of an iteration block
   specified by start iteration and end iteration (inclusive) to be iterated.
   The function then returns true if this chunk should be iterated
   and returns false if the end of the loop is reached (and the start/end
   iteration then contain useless values)
 */
typedef bool (scheduling_function_t) (int *start, int *end, int n, int p, int iam);


/*============================================================================*/
/* functions */

/* sequential scheduling */
void sequential_scheduling_setup (int n, int p, int iam);
bool sequential_scheduling (int *start_iteration, int *end_iteration, int n, int p, int iam);

/* block scheduling */
void block_scheduling_setup (int n, int p, int iam);
bool block_scheduling (int *start_iteration, int *end_iteration, int n, int p, int iam);

/* block-cyclic scheduling */
void blockcyclic_scheduling_setup (int n, int p, int iam);
bool blockcyclic_scheduling (int *start_iteration, int *end_iteration, int n, int p, int iam);

/* self scheduling */
void self_scheduling_setup (int n, int p, int iam);
bool self_scheduling (int *start_iteration, int *end_iteration, int n, int p, int iam);

/* guided self scheduling */
void gss_setup (int n, int p, int iam);
bool gss (int *start_iteration, int *end_iteration, int n, int p, int iam);

/* factoring */
void factoring_setup (int n, int p, int iam);
bool factoring (int *start_iteration, int *end_iteration, int n, int p, int iam);


#endif

/*============================================================================*
 *                             that's all folks                               *
 *============================================================================*/
