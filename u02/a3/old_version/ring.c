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
  int err, size, rank, res;
  MPI_Status status;
  // Init MPI
  err = MPI_Init( &argc, &argv);
  // get number of processors
  err = MPI_Comm_size(MPI_COMM_WORLD, &size);
  // get own number
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int a[size];
  int MASTER = 0;
  int LAST = size-1;

  if(rank == MASTER) {
    //printf("Master is sending now with rank: %d\n", rank);
    // first process
    a[rank] = rank+1; // add 1 at position 0
    // send to second process
    res = MPI_Send(a, // Sendepuffer (in)
                    size, // Anzahl Datenelemente (in)
                    MPI_INT, // MPI Datentyp (in)
                    1, // Empfänger (in)
                    1, // zusätzlicher Tag (in)
                    MPI_COMM_WORLD // Communicator (in)
    );
    if(res != MPI_SUCCESS){
      printf("Error sending from master to second process.\n");
      return 1;
    }
    // recieve filled vector from last process
    res = MPI_Recv(a, // Empfangspuffer (out)
                    size, // Anzahl Datenelemente (in)
                    MPI_INT, // MPI Datentyp (in)
                    LAST, // Sender (in) Last process
                    MPI_ANY_TAG, // zusätzlicher Tag (in)
                    MPI_COMM_WORLD, // Communicator (in)
                    &status // Empfangsstatus (out)
    );

    int CHECKSUM = createChecksum(size); // create checksum for result test
    
    int arraySum = createArraySum(a, size); // create sum of the filled vector
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
  } else {
    // middle processes
    //printf("Middle one is recieving and sending now with rank: %d", rank);
    res = MPI_Recv(a, // Empfangspuffer (out)
                    size, // Anzahl Datenelemente (in)
                    MPI_INT, // MPI Datentyp (in)
                    rank-1, // Sender (in) Last process
                    MPI_ANY_TAG, // zusätzlicher Tag (in)
                    MPI_COMM_WORLD, // Communicator (in)
                    &status // Empfangsstatus (out)
    );
    if(res != MPI_SUCCESS){
      printf("Error recieving from previous sender %d.\n", rank-1);
      return 1;
    }
    a[rank] = rank+1; // add rank number +1 to the array at my array position
    res = MPI_Send(a, // Sendepuffer (in)
                    size, // Anzahl Datenelemente (in)
                    MPI_INT, // MPI Datentyp (in)
                    ((rank+1)%size), // Empfänger (in)
                    1, // zusätzlicher Tag (in)
                    MPI_COMM_WORLD // Communicator (in)
    );
    if(res != MPI_SUCCESS){
      printf("Failed to send to next process %d.\n", rank+1);
      return 1;
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
