/*==============================================================================
  
   Purpose:    heat distribution
   Author:     Rudolf Berrendorf
               Computer Science Department
               Bonn-Rhein-Sieg University of Applied Sciences
	       53754 Sankt Augustin, Germany
               rudolf.berrendorf@h-brs.de

==============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include <libFHBRS.h>


/*============================================================================*/
// defines and global variables

#define DISPLAY_SIZE  1024	       // display window size
#define HEAT_MIN 20.0		       // minimum heat value
#define HEAT_MAX 100.0		       // maximum heat value

// display graphically?
static int graphics = 0;
static int window_number = -1;


/*============================================================================*/
// initialize h-array

static void
init_data (int xsize, int ysize, double h[xsize + 2][ysize + 2]) {
  // initialize data
  for (int i = 0; i < xsize + 2; i++)
    for (int j = 0; j < ysize + 2; j++)
      h[i][j] = HEAT_MIN;
}

/*============================================================================*/
// permanent heat source

static void
heat_source (int xsize, int ysize, double h[xsize + 2][ysize + 2]) {
  int sizex, sizey;

  // position 4 permanent heat sources
  sizex = xsize / 4;
  sizey = ysize / 4;

  // time-invariant constant values
  for (int i = 0; i < sizex; i++)
    for (int j = 0; j < sizey; j++) {
      h[1 * xsize / 4 - (sizex / 2) + i][1 * ysize / 4 - (sizey / 2) + j] = HEAT_MAX;
      h[1 * xsize / 4 - (sizex / 2) + i][3 * ysize / 4 - (sizey / 2) + j] = HEAT_MAX;
      h[3 * xsize / 4 - (sizex / 2) + i][1 * ysize / 4 - (sizey / 2) + j] = HEAT_MAX;
      h[3 * xsize / 4 - (sizex / 2) + i][3 * ysize / 4 - (sizey / 2) + j] = HEAT_MAX;
    }
}

/*============================================================================*/
// finite difference method (Laplace)

static void
update (int xsize, int ysize, double h[xsize + 2][ysize + 2], double h_new[xsize + 2][ysize + 2]) {
  // calculate new value based upon neighbor values
  for (int i = 1; i <= xsize; i++)
    for (int j = 1; j <= ysize; j++)
      h_new[i][j] = 0.25 * (h[i - 1][j] + h[i + 1][j] + h[i][j - 1] + h[i][j + 1]);

  // update array with new values
  for (int i = 1; i <= xsize; i++)
    for (int j = 1; j <= ysize; j++)
      h[i][j] = h_new[i][j];
}


/*============================================================================*/
// display current heat distribution

static void
display (int xsize, int ysize, double h[xsize + 2][ysize + 2]) {
  for (int i = 1; i <= xsize; i++)
    for (int j = 1; j <= ysize; j++) {
      double color;
      
      color = (h[i][j] - HEAT_MIN) / (HEAT_MAX - HEAT_MIN);
      
      if (color < 0.0)
        color = 0.0;
      if (color > 1.0)
        color = 1.0;
      
      if (graphics) {
        graphic_setRainbowColor (window_number, (double) color);
        graphic_drawPoint (window_number, i, j);
      }
    }
  
  // display columnwise
  if (graphics)
    graphic_flush (window_number);
}


/*============================================================================*/
// calculate checksum

static unsigned long
checksum (int xsize, int ysize, double h[xsize + 2][ysize + 2]) {
  unsigned long checksum = 0;

  for (int i = 0; i < xsize + 2; i++)
    for (int j = 0; j < ysize + 2; j++)
      checksum += (unsigned long) h[i][j];

  return checksum;
}

/*============================================================================*/

static void
usage (char *progname) {
  printf ("usage: %s size niter graphics\n", progname);
  printf ("where\tsize mesh size in each dimension (e.g. 2000)\n");
  printf ("\tniter number of iterations (e.g. 1000)\n");
  printf ("\tdisplay graphically every x steps (e.g. 100). 0 means no graphic\n");
  exit (EXIT_FAILURE);
}


/*============================================================================*/

int
main (int argc, char **argv) {
  int niter, size;
  double t0, t1, ttotal;


  /*--------------------------------------------------------------------------*/
  // check optional runtime parameters

  if (argc == 4) {
    if (sscanf (argv[1], "%d", &size) != 1)
      usage (argv[0]);
    if (sscanf (argv[2], "%d", &niter) != 1)
      usage (argv[0]);
    if (sscanf (argv[3], "%d", &graphics) != 1)
      usage (argv[0]);
  } else {
    // error
    usage (argv[0]);
  }

  /*--------------------------------------------------------------------------*/

  // start graphics
  if (graphics) {
    if ((window_number = graphic_start (DISPLAY_SIZE, DISPLAY_SIZE, "Heat Distribution")) < 0)
      return EXIT_FAILURE;
    graphic_userCoordinateSystem (window_number, 0.0, 0.0, size, size);
  }

  // finite difference array with additional 2 boundary cells to simplify calculation
  double *h_buf = malloc ((size + 2) * (size + 2) * sizeof (*h_buf));
  double *h_new_buf = malloc ((size + 2) * (size + 2) * sizeof (*h_new_buf));
  if ((h_buf == NULL) || (h_new_buf == NULL)) {
    // error
    printf ("no more memory\n");
    exit (EXIT_FAILURE);
  }
  
  int xsize = size;
  int ysize = size;

  // for better reading/casting define an array type and appropriate variables
  typedef double array_t[ysize + 2];
  array_t *h = (array_t *) h_buf;
  array_t *h_new = (array_t *) h_new_buf;

  // initialize data
  init_data (xsize, ysize, h);

  // add heat sources
  heat_source (xsize, ysize, h);

  /*--------------------------------------------------------------------------*/

  // iteration steps
  t0 = ttotal = gettime ();
  if (graphics)
    display (xsize, ysize, h);

  for (int t = 0; t < niter; t++) {
    // display every niter step
    if ((graphics > 0) && (t % graphics == (graphics - 1))) {
      t1 = gettime ();
      display (xsize, ysize, h);
      printf ("time step %8d of %8d: %8.6f s per time step\n", (t + 1), niter, (t1 - t0) / graphics);
      t0 = gettime ();
    }

    // compute new values
    update (xsize, ysize, h, h_new);
    
    // permanent heat source
    heat_source (xsize, ysize, h);
  }

  /*--------------------------------------------------------------------------*/

  ttotal = gettime () - ttotal;
  printf ("total time: %.6f\n", ttotal);
  printf ("checksum  : %lu\n", checksum (xsize, ysize, h));

  // finish graphics
  if (graphics)
    graphic_end (window_number);

  free (h);
  free (h_new);

  return EXIT_SUCCESS;
}


/*============================================================================*
 *                             that's all folks                               *
 *============================================================================*/
