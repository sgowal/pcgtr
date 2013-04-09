/*****************************************
 * Optimization of the sphere function:  *
 *   min f(x) = \sum_{i = 1}^N x_i^2     *
 * We removed all the unnecessary bits!  *
 *****************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Include PCG solver.
#include <pcgtr.h>

// Declaration of the cost function.
float MyCost(matrix_t *x, matrix_t *J, matrix_t *H, void *p);

int main(int argc, char const *argv[]) {
  float *values;
  int i, n;
  pcgtr_solver_t solver;

  // Read number of variables to optimize from commandline.
  n = 10;
  if (argc > 1) {
    n = atoi(argv[1]);
    if (n <= 1) {
      fprintf(stderr, "\nUsage: %s [N]\n"
                      "  where N is the number of variable to optimize (N > 1).\n\n", argv[0]);
      return -1;
    }
  }
  // Initialize storage.
  values = malloc(n*sizeof(float));
  // Create initial estimate.
  for (i = 0; i < n; ++i)
    values[i] = rand()/(float)RAND_MAX;
  matrix_t x = {n, 1, values};

  // Initialize solver.
  pcgtr_solver_init(&solver, n, 0);

  // Set functions -- after the initialization.
  solver.functions.f = MyCost;
  // Set settings.
  solver.settings.verbose = 1;

  // Solve.
  pcgtr_solver_solve(&solver, &x);
  
  // Destroy. 
  pcgtr_solver_destroy(&solver);

  // Clear storage.
  free(values);

  return 0;
}

float MyCost(matrix_t *x, matrix_t *J, matrix_t *H, void *p) {
  int i;
  float c = 0.0f;

  // Weighted sphere function.
  for (i = 0; i < x->m; ++i) {
    float w = (float)i;
    c += w*x->v[i]*x->v[i];
    J->v[i] = 2.0f*w*x->v[i];
    H->v[matrix_index(H,i,i)] = 2.0f*w; 
  }

  return c;
}