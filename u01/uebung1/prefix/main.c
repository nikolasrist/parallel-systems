//===================================================================================
// apply a prefix operation

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <libFHBRS.h>

#include "prefix.h"


//===================================================================================
// prefix operator instance (long addition)

static void long_add(void* r, void* a, void* b)
{
  *((long*) r) = *((long*) a) + *((long*) b);
}


//===================================================================================

int main(int argc, char **argv)
{
  prefix_data data;
  prefix_operator op;
  int count;
  // number of measurements
  int n_times = (argc == 1) ? 3 : atoi(argv[1]);
  double t[n_times];
  double t_min = 1e38;


  // allocate vector
  long n = 1 << 29;
  long *buffer = malloc(n * sizeof(*buffer));
  assert(buffer != NULL);
  
  // initialize vector data structure
  data.pointer = (char *)buffer;
  data.bytes_len = sizeof(buffer[0]);
  data.vector_len = n;
  
  // initialize operator structure
  int null_value = 0;
  op.name = "add";
  op.function = long_add;
  op.neutral = &null_value;
  op.type_len = sizeof(buffer[0]);

  
  // repeat measurements
  (void) gettime();
  for(int iter=0; iter<n_times; iter++)
    {
      // initialize argument vector
      for(long i=0; i<n; i++)
	buffer[i] = i+1;

      // call prefix operation
      t[iter] = gettime();
      int rc = prefix(data, op, &count);
      t[iter] = gettime() - t[iter];
      assert(rc == PREFIX_SUCCESS);

      // check for minimum time
      if(t[iter] < t_min)
	t_min = t[iter];

      // check correctness of result
      long sum = 0;
      for(long i=0; i<n; i++)
	{
	  sum += (i+1);
	  if(buffer[i] != sum)
	    {
	      printf("first error at position %ld (%ld should be %ld)\n", i, buffer[i], sum);
	      break;
	    }
	}
    }


  // print minimum time consumed for one prefix operation
  printf("minimum time for prefix operation: %.6f s\n", t_min);


  // clean up
  free(buffer);

  return 0;
}

//===================================================================================
