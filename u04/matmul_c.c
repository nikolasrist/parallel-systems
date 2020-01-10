/*==============================================================================
  
   Purpose          : matrix multiply
   Author           : Rudolf Berrendorf
                      Computer Science Department
                      Bonn-Rhein-Sieg University of Applied Sciences
	              53754 Sankt Augustin, Germany
                      rudolf.berrendorf@h-brs.de
  
==============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(MKL)
#include <mkl.h>
#else
#include <cblas.h>
#define dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) \
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, *m, *n, *k, *alpha, a, *n, b, *m, *beta, c, *k)
#endif

#include <libFHBRS.h>

/*===========================================================================*/

// switch between production version (large matrix) and test version (small matrix)
#if 1
#define N 2999            // size of matrix
#else
#define N 5
#endif

#define EPS 1e-4          // tolerance on correctness proof


// three matrix sets
static double a1[N][N], b1[N][N], c1[N][N];
static double a2[N][N], b2[N][N], c2[N][N];
static double a3[N][N], b3[N][N], c3[N][N];

// forward declaration
static int matprint(char *header, int n, double a[n][n]);

/*===========================================================================*/
// apply your optimizations to this version

static void matmul(int n, double a[n][n], double b[n][n], double c[n][n]) {
  int t = 272;
  double tmp;
  double j21, j22, j23, j24, j25, j26, j27, j28, j29, j210, j211, j212, j213, j214, j215;
  int i, i2, j, j2, k, k2, iUpper, kUpper, jUpper, iBound, kBound, jBound, jBound2;
  for (i = 0; i < n; i += t){
    iUpper=i+t;
    iBound = iUpper < n ? iUpper : n;
    for (k = 0; k < n; k += t){
      kUpper=k+t;
      kBound = kUpper < n ? kUpper : n;
      for (j = 0; j < n; j += t){
        jUpper=j+t;
        jBound = jUpper < n ? jUpper : -1;
        jBound2 = jUpper > n ? n : -1;
        for (i2 = i; i2 < iBound; ++i2)
          for (k2 = k; k2 < kBound; ++k2){
            tmp = b[i2][k2];
            for (j2 = j; j2 < jBound;j2+=16){
              a[i2][j2] +=  tmp * c[k2][j2];
              a[i2][j2+1] +=  tmp * c[k2][j2+1];
              a[i2][j2+2] +=  tmp * c[k2][j2+2];
              a[i2][j2+3] +=  tmp * c[k2][j2+3];
              a[i2][j2+4] +=  tmp * c[k2][j2+4];
              a[i2][j2+5] +=  tmp * c[k2][j2+5];
              a[i2][j2+6] +=  tmp * c[k2][j2+6];
              a[i2][j2+7] +=  tmp * c[k2][j2+7];
              a[i2][j2+8] +=  tmp * c[k2][j2+8];
              a[i2][j2+9] +=  tmp * c[k2][j2+9];
              a[i2][j2+10] +=  tmp * c[k2][j2+10];
              a[i2][j2+11] +=  tmp * c[k2][j2+11];
              a[i2][j2+12] +=  tmp * c[k2][j2+12];
              a[i2][j2+13] +=  tmp * c[k2][j2+13];
              a[i2][j2+14] +=  tmp * c[k2][j2+14];
              a[i2][j2+15] +=  tmp * c[k2][j2+15];
            }
            for (j2 = j; j2 < jBound2; ++j2){
              a[i2][j2] +=  tmp * c[k2][j2];
            }
          }
      }
    }
  }
}

/*===========================================================================*/
// original version

static void matmul_original(int n, double a[n][n], double b[n][n], double c[n][n]) {
  for(int j=0; j<n; j++)
    for(int k=0; k<n; k++)
      for(int i=0; i<n; i++)
	a[i][j] += b[i][k] * c[k][j];
}


/*===========================================================================*/
// highly optimized version

static void matmul_optimized(int n, double a[n][n], double b[n][n], double c[n][n]) {
  const char trans = 'N';
  const double alpha = 1.0;
  const double beta = 0.0;

  // call BLAS DGEMM ( C interface)
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              n, n, n, alpha, (double *)b, n, (double *)c, n, beta, (double *)a, n);
}


/*===========================================================================*/
// initialize matrices

static void matinit(int n, double a[n][n], double b[n][n], double c[n][n]) {
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++) {
      a[i][j] = 0.0;
      b[i][j] = 2*j - i;
      c[i][j] = 5*i + j;
    }
}


/*===========================================================================*/
// check for correctness against original version

static int matcheck(int n, double a1[n][n], double a2[n][n]) {
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      if(fabs(a1[i][j] - a2[i][j]) > EPS) {
        printf("error at position a[%d][%d]: %f should be %f\n", i, j, a1[i][j], a2[i][j]);
        return -1;
      }
     
  return 0;
}


/*===========================================================================*/
// print small matrices to stdout

static int matprint(char *header, int n, double a[n][n]) {
  printf("%s\n", header);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      printf("%4.1f ", a[i][j]);
    }
    printf("\n");
  }
}

/*===========================================================================*/

int main(int argc, char **argv) {

  double t1, t2, t3;
  int n = N;


  // start optimized version
  matinit(n, a1, b1, c1);
  t1 = gettime();
  matmul_optimized(n, a1, b1, c1);
  t1 = gettime() - t1;
  printf("%10s: %9.3f s, %9.3f GFLOPS\n", "optimized", t1, 2.0*n*n*n/1e9/t1);


  // start your version
  matinit(n, a2, b2, c2);
  t2 = gettime();
  matmul(n, a2, b2, c2);
  t2 = gettime() - t2;
  // check results for correctness
  if(matcheck(n, a1, a2) == 0) {
    printf("%10s: %9.3f s, %9.3f GFLOPS\n", "your", t2, 2.0*n*n*n/1e9/t2);
  } else {
    printf("incorrect result for your version\n");
    if(n < 10) {
      matprint("matrix b", n, b2);
      matprint("matrix c", n, c2);
      matprint("matrix a", n, a2);
      matprint("matrix should be", n, a1);
    }
  }


  // start original version
  matinit(n, a3, b3, c3);
  t3 = gettime();
  matmul_original(n, a3, b3, c3);
  t3 = gettime() - t3;
  // check results for correctness
  if(matcheck(n, a1, a3) == 0) {
    printf("%10s: %9.3f s, %9.3f GFLOPS\n", "original", t3, 2.0*n*n*n/1e9/t3);
  } else {
    printf("incorrect result for optimized version\n");
    if(n < 10) {
      matprint("matrix b", n, b2);
      matprint("matrix c", n, c2);
      matprint("matrix a", n, a2);
      matprint("matrix should be", n, a1);
    }
  }

  return 0;
}

/*============================================================================*
 *                             that's all folks                               *
 *============================================================================*/
