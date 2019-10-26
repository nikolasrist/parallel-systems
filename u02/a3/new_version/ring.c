/*==============================================================================
  
   Purpose          : MPI pattern
   Author           : Rudolf Berrendorf
                      Computer Science Department
                      Bonn-Rhein-Sieg University of Applied Sciences
	              53754 Sankt Augustin, Germany
                      rudolf.berrendorf@h-brs.de
  
==============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <libFHBRS.h>
#include <mpi.h>

//==============================================================================

int createChecksum(int n);
int createArraySum(int a[], int size);

int main(int argc, char **argv)
{ 
  // create initial array to send to next process
  int err, size, rank;
  MPI_Status status;
  // Init MPI
  err = MPI_Init( &argc, &argv);
  // get number of processors
  err = MPI_Comm_size(MPI_COMM_WORLD, &size);
  // get own number
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int MASTER = 0;
  void *gather_array = NULL;
  if(rank == MASTER) {
    gather_array = (int *)malloc(sizeof(int) * size);
  }

  // increase rank by one to start from 1
  int array_value = rank+1;

  MPI_Gather(&array_value, 1, MPI_INT, gather_array,1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if(rank == MASTER) {
    int CHECKSUM = createChecksum(size); // create checksum for result test
    

    int arraySum = createArraySum(gather_array, size); // create sum of the filled vector
    if(CHECKSUM == arraySum) { 
      // each process added correctly his rank number +1
      printf("Successful run.\n");
      printf("Checksum: %d\n", CHECKSUM);
      printf("ArraySum: %d\n", arraySum);
    } else {
      // any process added a wrong value
      printf("Failed to run program.\n");
      printf("Checksum: %d\n", CHECKSUM);
      printf("ArraySum: %d\n", arraySum);
    }
  }
  err = MPI_Finalize();
  return 0;
}

int createArraySum (int a[], int size) {
  int result = 0;
  printf("Result-Vector:\n");
  for (int i=0; i < size; i++) {
    printf("%d,", a[i]);
    result += a[i];
  }
  return result;
}

int createChecksum(int n) {
  int result = 0;
  for (int i = 1; i <= n; i++){
    result += i;
  }
  return result;
}

/*============================================================================*
 *                             that's all folks                               *
 *============================================================================*/
