/*************************************************************************
 * Optimization of the Rosenbrock function:                              *
 *   min f(x) = \sum_{i = 1}^{N-1} 100 (x_{i+1} - x_i^2)^2 + (x_i - 1)^2 *
 *************************************************************************/

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
  float *values, *correct;
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
  // Solution of the Rosenbrock function.
  correct = malloc(n*sizeof(float));
  for (i = 0; i < n; ++i)
    correct[i] = 1.0f;
  matrix_t correct_x = {n, 1, correct};

  // Initialize solver.
  pcgtr_solver_init(&solver, n, sparse);

  // Set functions.
  solver.functions.f = MyCost;
  if (sparse) {
    solver.functions.HessianMultiplication = MyHessianMultiplication;
    solver.functions.DiagonalInverseHessianApproximation = MyDiagonalInverseHessianApproximation;
    pcgtr_solver_halloc(&solver, n, 3);
  }

  // Set settings.
  solver.settings.verbose = (trials == 1);

  total_time = 0.0f;
  total_cost = 0.0f;
  total_distance = 0.0f;
  for (k = 0; k < trials; ++k) {
    // Create initial estimate around the solution.
    for (i = 0; i < n; ++i)
      values[i] = 2.0f*(rand()/(float)RAND_MAX);
    matrix_t x = {n, 1, values};
    if (trials == 1)
      matrix_print(&x, "xinit");

    // Solve.
    current_time = gettime();
    pcgtr_solver_solve(&solver, &x);
    total_time += gettime() - current_time;
    total_cost += MyCost(&x, NULL, NULL, NULL);
    if (trials == 1)
      matrix_print(&x, "xfinal");
    else
      printf(".");
    fflush(stdout);
    total_distance += vector_norm(matrix_sub(&x, &correct_x, &x)); // This call modifies x.
  }
  free(values);
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
  int n = x->m;

  for (i = 0; i < x->m-1; ++i)
    c += 100.0f * (x->v[i+1] - x->v[i]*x->v[i]) * (x->v[i+1] - x->v[i]*x->v[i]) + (x->v[i] - 1.0f) * (x->v[i] - 1.0f);

  // Avoid if-statement within for-loop.
  if (J) {
    J->v[0] = - 400.0f * x->v[0] * (x->v[1] - x->v[0]*x->v[0]) + 2.0f * (x->v[0] - 1.0f);
    for (i = 1; i < x->m-1; ++i)
      J->v[i] = 200.0f * (x->v[i] - x->v[i-1]*x->v[i-1]) - 400.0f * x->v[i] * (x->v[i+1] - x->v[i]*x->v[i]) + 2.0f * (x->v[i] - 1.0f);
    J->v[n-1] = 200.0f * (x->v[n-1] - x->v[n-2]*x->v[n-2]);
  }

  // Compute Hessian depending on whether in sparse mode or not.
  // Avoid if-statement within for-loop.
  if (H) {
    if (H->n == H->m) {
      H->v[matrix_index(H,0,0)] = 2.0f + 1200.0f * x->v[0] * x->v[0] - 400.0f * x->v[1];
      H->v[matrix_index(H,0,1)] = -400.0f * x->v[0];
      for (i = 1; i < n-1; ++i) {
        H->v[matrix_index(H,i,i-1)] = -400.0f * x->v[i-1];
        H->v[matrix_index(H,i,i)] = 202.0f + 1200.0f * x->v[i] * x->v[i] - 400.0f * x->v[i+1];
        H->v[matrix_index(H,i,i+1)] = -400.0f * x->v[i];
      }
      H->v[matrix_index(H,n-1,n-2)] = - 400.0f * x->v[n-2];
      H->v[matrix_index(H,n-1,n-1)] = 200.0f;
    } else {
      H->v[matrix_index(H,0,0)] = 2.0f + 1200.0f * x->v[0] * x->v[0] - 400.0f * x->v[1];
      H->v[matrix_index(H,0,1)] = -400.0f * x->v[0];
      for (i = 1; i < n-1; ++i) {
        H->v[matrix_index(H,i,0)] = -400.0f * x->v[i-1];
        H->v[matrix_index(H,i,1)] = 202.0f + 1200.0f * x->v[i] * x->v[i] - 400.0f * x->v[i+1];
        H->v[matrix_index(H,i,2)] = -400.0f * x->v[i];
      }
      H->v[matrix_index(H,n-1,1)] = - 400.0f * x->v[n-2];
      H->v[matrix_index(H,n-1,2)] = 200.0f;
    }
  }

  return c;
}

#define MAX(x,y) (((x)>(y))?(x):(y))
matrix_t *MyDiagonalInverseHessianApproximation(matrix_t *a, matrix_t *z, void *p) {
  int i;
  int m = a->m;

  for (i = 0; i < m; i++)
    z->v[i] = 1.0f/MAX(sqrtf(a->v[matrix_index(a,i,0)]*a->v[matrix_index(a,i,0)] +
                             a->v[matrix_index(a,i,1)]*a->v[matrix_index(a,i,1)] +
                             a->v[matrix_index(a,i,2)]*a->v[matrix_index(a,i,2)]), 0.01f);

  return z;
}

matrix_t *MyHessianMultiplication(matrix_t *a, matrix_t *b, matrix_t *z, void *p) {
  int i, j;
  int m = a->m;
  int o = b->n;

  for (j = 0; j < o; j++) {
    z->v[matrix_index(z,0,j)] = a->v[matrix_index(a,0,0)]*b->v[matrix_index(b,0,j)] +
                                a->v[matrix_index(a,0,1)]*b->v[matrix_index(b,1,j)];
    for (i = 1; i < m-1; i++) {
      z->v[matrix_index(z,i,j)] = a->v[matrix_index(a,i,0)]*b->v[matrix_index(b,i-1,j)] +
                                  a->v[matrix_index(a,i,1)]*b->v[matrix_index(b,i,j)] +
                                  a->v[matrix_index(a,i,2)]*b->v[matrix_index(b,i+1,j)];
    }
    z->v[matrix_index(z,m-1,j)] = a->v[matrix_index(a,m-1,1)]*b->v[matrix_index(b,m-2,j)] +
                                  a->v[matrix_index(a,m-1,2)]*b->v[matrix_index(b,m-1,j)];
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