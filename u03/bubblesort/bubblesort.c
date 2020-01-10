/*==============================================================================

   Purpose          : bubblesort
   Author           : Rudolf Berrendorf
                      Computer Science Department
                      Bonn-Rhein-Sieg University of Applied Sciences
                      53754 Sankt Augustin, Germany
                      rudolf.berrendorf@h-brs.de

==============================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <libFHBRS.h>

// type for data
typedef double data_t;

//==============================================================================

// swap content of two variables
#define swap(x,y) { data_t tmp = *x; *x=*y; *y=tmp;}

//==============================================================================

static void bubblesort(int n, data_t a[n])
{
  for (int i = 1; i < n; i++)
    for (int j = 0; j < n - i; j++)
      if (a[j] > a[j + 1])
        swap (&a[j], &a[j + 1]);

}

//==============================================================================

int main(int argc, char **argv)
{
  int n = atoi(argv[1]);
  data_t *a = malloc(n * sizeof(*a));
  if(a == NULL)
    {
      printf("no more memory\n");
      exit(1);
    }


  // initialize data
  srand(0);
  for(int i=0; i<n; i++)
    a[i] = rand();

  // sort
  bubblesort(n, a);

  // check result
  for(int i=0; i<n-1; i++)
    assert(a[i] <= a[i+1]);

  return 0;
}

/*============================================================================*
 *                             that's all folks                               *
 *============================================================================*/
