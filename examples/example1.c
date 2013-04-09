/*****************************************
 * Optimization of the sphere function:  *
 *   min f(x) = \sum_{i = 1}^N x_i^2     *
 *****************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

// Include PCG solver.
#include <pcgtr.h>

// Cost function declaration.
// This function needs to return the cost at x and should
// compute the Jacobian J and Hessian H.
// It may take some parameters through p.
float MyCost(matrix_t *x, matrix_t *J, matrix_t *H, void *p);
matrix_t *MyHessianMultiplication(matrix_t *a, matrix_t *b, matrix_t *z, void *p);
matrix_t *MyDiagonalInverseHessianApproximation(matrix_t *a, matrix_t *z, void *p);
// Helper functions
double gettime();
void seedrand();

int main(int argc, char const *argv[]) {
  float *values, *weights;
  double total_time, current_time, total_cost, total_distance;
  int i, k, n, trials, sparse;
  pcgtr_solver_t solver;

  // Initialize random number generator (not needed).
  seedrand();

  // Read number of variables to optimize from commandline.
  n = 10;
  trials = 1;
  sparse = 0;
  if (argc > 1) {
    n = atoi(argv[1]);
    if (argc > 2)
      trials = atoi(argv[2]);
    if (argc > 3)
      sparse = (strncmp(argv[3], "sparse", 6) == 0);
    if (n <= 1 || trials <= 0) {
      fprintf(stderr, "\nUsage: %s [N] [trials] [sparse]\n"
                      "  where N is the number of variable to optimize (N > 1),\n"
                      "        trials is the number of trials to perform (Trials > 0),\n"
                      "        sparse is set to \"sparse\" to indicate sparse processing.\n\n", argv[0]);
      return -1;
    }
  }
  // Initialize storage.
  values = malloc(n*sizeof(float));
  weights = malloc(n*sizeof(float));

  // Initialize solver.
  pcgtr_solver_init(&solver, n, sparse);

  // Set functions.
  solver.functions.f = MyCost;
  if (sparse) {
    solver.functions.HessianMultiplication = MyHessianMultiplication;
    solver.functions.DiagonalInverseHessianApproximation = MyDiagonalInverseHessianApproximation;
    pcgtr_solver_halloc(&solver, n, 1);
  }
  solver.functions.p = weights;

  // Set settings.
  solver.settings.verbose = (trials == 1);

  total_time = 0.0f;
  total_cost = 0.0f;
  total_distance = 0.0f;
  for (k = 0; k < trials; ++k) {
    // Create initial estimate.
    for (i = 0; i < n; ++i)
      values[i] = rand()/(float)RAND_MAX;
    matrix_t x = {n, 1, values};
    if (trials == 1)
      matrix_print(&x, "xinit");

    // Create a random weighting for the cost function.
    for (i = 0; i < n; ++i)
      weights[i] = rand()/(float)RAND_MAX + 0.1f;

    // Solve.
    current_time = gettime();
    pcgtr_solver_solve(&solver, &x);
    total_time += gettime() - current_time;
    total_cost += MyCost(&x, NULL, NULL, weights);
    total_distance += vector_norm(&x);
    if (trials == 1)
      matrix_print(&x, "xfinal");
    else
      printf(".");
    fflush(stdout);
  }
  free(values);
  free(weights);
  printf("\nAverage time taken: %.2f ms (%0.2f us)\n", 1000.0f*total_time/trials, 1000000.0f*total_time/trials);
  printf("Average distance to correct solution: %g\n", total_distance/trials);
  printf("Average final cost: %g\n\n", total_cost/trials);

  // Destroy. 
  pcgtr_solver_destroy(&solver);

  return 0;
}

float MyCost(matrix_t *x, matrix_t *J, matrix_t *H, void *p) {
  int i;
  float c = 0.0f;
  float *w = (float *)p;


  for (i = 0; i < x->m; ++i)
    c += w[i]*x->v[i]*x->v[i];

  // Avoid if-statement within for-loop.
  if (J) {
    for (i = 0; i < x->m; ++i)
      J->v[i] = 2.0f*w[i]*x->v[i];
  }

  // Compute Hessian depending on whether in sparse mode or not.
  // Avoid if-statement within for-loop.
  if (H) {
    if (H->n == H->m) {
      for (i = 0; i < x->m; ++i)
        H->v[matrix_index(H,i,i)] = 2.0f*w[i];
    } else {
      for (i = 0; i < x->m; ++i)
        H->v[i] = 2.0f*w[i];
    }
  }

  return c;
}

matrix_t *MyDiagonalInverseHessianApproximation(matrix_t *a, matrix_t *z, void *p) {
  int i;
  int m = a->m;

  for (i = 0; i < m; i++) {
    z->v[i] = 1.0f/a->v[i];
  }

  return z;
}

matrix_t *MyHessianMultiplication(matrix_t *a, matrix_t *b, matrix_t *z, void *p) {
  int i, j;
  int m = a->m;
  int o = b->n;

  for (i = 0; i < m; i++) {
    for (j = 0; j < o; j++) {
      z->v[i*o+j] = a->v[i]*b->v[i*o+j];
    }
  }

  return z;
}

double gettime() {
  struct timeval time;
  gettimeofday(&time, NULL);
  return (double)time.tv_sec + ((double)time.tv_usec)/1000000.0;
}

void seedrand() {
  struct timeval time;
  gettimeofday(&time, NULL);
  srand(time.tv_usec);
}