/*==============================================================================
  
   Purpose          : matrix multiply
   Author           : Rudolf Berrendorf
                      Department of Computer Science
                      Bonn-Rhein-Sieg University of Applied Sciences
	              53754 Sankt Augustin, Germany
                      rudolf.berrendorf@h-brs.de
  
==============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <libFHBRS.h>

//==============================================================================

static void matmul_ikj(int n, double a[n][n], double b[n][n], double c[n][n])
{
  for(int i=0; i<n; i++)
    for(int k=0; k<n; k++)
      for(int j=0; j<n; j++)
	{
	  a[i][j] += b[i][k] * c[k][j];
	}
}

//==============================================================================

static void matmul_ijk(int n, double a[n][n], double b[n][n], double c[n][n])
{
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      for(int k=0; k<n; k++)
	{
	  a[i][j] += b[i][k] * c[k][j];
	}
}

//==============================================================================

static void matmul_jik(int n, double a[n][n], double b[n][n], double c[n][n])
{
  for(int j=0; j<n; j++)
    for(int i=0; i<n; i++)
      for(int k=0; k<n; k++)
	{
	  a[i][j] += b[i][k] * c[k][j];
	}
}

//==============================================================================

static void matmul_jki(int n, double a[n][n], double b[n][n], double c[n][n])
{
  for(int j=0; j<n; j++)
    for(int k=0; k<n; k++)
      for(int i=0; i<n; i++)
	{
	  a[i][j] += b[i][k] * c[k][j];
	}
}

//==============================================================================

static void matmul_kij(int n, double a[n][n], double b[n][n], double c[n][n])
{
  for(int k=0; k<n; k++)
    for(int j=0; j<n; j++)
      for(int i=0; i<n; i++)
	{
	  a[i][j] += b[i][k] * c[k][j];
	}
}

//==============================================================================

static void matmul_kji(int n, double a[n][n], double b[n][n], double c[n][n])
{
  for(int k=0; k<n; k++)
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
	{
	  a[i][j] += b[i][k] * c[k][j];
	}
}

//==============================================================================

int main(int argc, char **argv)
{
  if(argc != 2)
    {
      printf("usage: %s n\n", argv[0]);
      exit(1);
    }

  int n = atoi(argv[1]);
  double a[n][n], b[n][n], c[n][n];
  
  //----------------------------------------------------------------------------


  //matmul_ijk(n, a, b, c);


  //----------------------------------------------------------------------------

  //matmul_ikj(n, a, b, c);

  //----------------------------------------------------------------------------

  matmul_jik(n, a, b, c);

  //----------------------------------------------------------------------------

  //matmul_jki(n, a, b, c);

  //----------------------------------------------------------------------------

  //matmul_kij(n, a, b, c);

  //----------------------------------------------------------------------------

  //matmul_kji(n, a, b, c);

  //----------------------------------------------------------------------------

  return 0;
}

/*============================================================================*
 *                             that's all folks                               *
 *============================================================================*/
