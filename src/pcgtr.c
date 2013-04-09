#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "pcgtr.h"

#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))

/********************
 * Helper functions *
 ********************/

static float step(pcgtr_solver_t *solver, matrix_t *x, float delta, int computeZ, int *positiveDefinite);
static void precondition(pcgtr_solver_t *solver, int *positiveDefinite);
static void trustregion(pcgtr_solver_t *solver, float delta);
static int get_solution(pcgtr_solver_t *solver, float lambda, float *eigenvalues, float *eigenvectors, float *alpha, float *s, int failsafe);
static float secular_equation(pcgtr_solver_t *solver, float lambda, float *eigenvalues, float *alpha, float delta);
static float rfzero(pcgtr_solver_t *solver, float x, float *eigenvalues, float *alpha, float delta, float *value);
static matrix_t *DefaultHessianMultiplication(matrix_t *a, matrix_t *b, matrix_t *z, void *p);
static matrix_t *DefaultDiagonalInverseHessianApproximation(matrix_t *a, matrix_t *z, void *p);

/**********************
 * The main functions *
 **********************/

int pcgtr_solver_init(pcgtr_solver_t *solver, int n, int sparse) {
  if (n < 2) {
    fprintf(stderr, "To avoid additional checks, cgsolve does not optimize less than 2 values.");
    return -1;
  }

  solver->settings.MaxPCGIter = n;
  solver->settings.MaxIter = 50;
  solver->settings.MaxZeroIter = 50;
  solver->settings.TolFun = 1e-6f;
  solver->settings.TolX = 1e-6f;
  solver->settings.TolZero = 1e-6f;
  solver->settings.eps = 1e-5f;
  solver->settings.verbose = 0;

  solver->functions.f = NULL;
  solver->functions.HessianMultiplication = DefaultHessianMultiplication;
  solver->functions.DiagonalInverseHessianApproximation = DefaultDiagonalInverseHessianApproximation;
  solver->functions.p = NULL;

  assert(matrix_init(&solver->grad, n, 1));
  assert(matrix_init(&solver->xOld, n, 1));
  assert(matrix_init(&solver->gradOld, n, 1));
  if (!sparse) {
    assert(matrix_init(&solver->hessian, n, n));
    assert(matrix_init(&solver->hessianOld, n, n));
  } else {
    solver->hessian.v = NULL;
    solver->hessianOld.v = NULL;
  }
  assert(matrix_init(&solver->Z, n, 2));
  assert(matrix_init(&solver->v1, n, 1));
  assert(matrix_init(&solver->v2, n, 1));
  assert(matrix_init(&solver->dx, n, 1));

  assert(matrix_init(&solver->trustG, 2, 1));
  assert(matrix_init(&solver->trustH, 2, 2));
  assert(matrix_init(&solver->A, n, 2));
  assert(matrix_init(&solver->trustS, 2, 1));
  assert(matrix_init(&solver->B, 2, 1));

  assert(matrix_init(&solver->rpcg, n, 1));
  assert(matrix_init(&solver->zpcg, n, 1));
  assert(matrix_init(&solver->ppcg, n, 1));
  assert(matrix_init(&solver->Appcg, n, 1));
  assert(matrix_init(&solver->Mpcg, n, 1));

  assert(matrix_init(&solver->eigenvaluesTR, 2, 1));
  assert(matrix_init(&solver->eigenvectorsTR, 2, 2));

  return 0;
}

void pcgtr_solver_halloc(pcgtr_solver_t *solver, int row, int col) {
  assert(matrix_init(&solver->hessian, row, col));
  assert(matrix_init(&solver->hessianOld, row, col));
}

void pcgtr_solver_destroy(pcgtr_solver_t *solver) {
  matrix_destroy(&solver->grad);
  matrix_destroy(&solver->xOld);
  matrix_destroy(&solver->gradOld);
  if (solver->hessian.v)
    matrix_destroy(&solver->hessian);
  if (solver->hessianOld.v)
    matrix_destroy(&solver->hessianOld);
  matrix_destroy(&solver->Z);
  matrix_destroy(&solver->v1);
  matrix_destroy(&solver->v2);
  matrix_destroy(&solver->dx);

  matrix_destroy(&solver->trustG);
  matrix_destroy(&solver->trustH);
  matrix_destroy(&solver->A);
  matrix_destroy(&solver->trustS);
  matrix_destroy(&solver->B);

  matrix_destroy(&solver->rpcg);
  matrix_destroy(&solver->zpcg);
  matrix_destroy(&solver->ppcg);
  matrix_destroy(&solver->Appcg);
  matrix_destroy(&solver->Mpcg);

  matrix_destroy(&solver->eigenvaluesTR);
  matrix_destroy(&solver->eigenvectorsTR);
}

void pcgtr_solver_solve(pcgtr_solver_t *solver, matrix_t *x) {
  int iter, n, positiveDefinite, done, computeZ;
  float f, normDx, delta, ratio, delbnd, fOld, normGrad, df, value;

  // Sanity checks.
  assert(x && x->v && (x->m == solver->grad.m));
  assert(solver->functions.f);
  assert(solver->functions.HessianMultiplication);
  assert(solver->functions.DiagonalInverseHessianApproximation);
  assert(solver->hessian.v && solver->hessianOld.v);

  n = x->m;

  // Initialization
  f = solver->functions.f(x, &solver->grad, &solver->hessian, solver->functions.p);
  if (solver->settings.verbose) {
    printf("Iteration        f(x)\n");
    printf("#%d           %10f\n", 0, f);
  }

  done = 0;
  fOld = 1.0f/0.0f; // infinity
  positiveDefinite = 1;
  normDx = 1.0f;
  delta = 10.0f;
  ratio = 0.0f;
  delbnd = MAX(100.0f*vector_norm(x), 1.0f);
  computeZ = 1;

  for (iter = 1; iter <= solver->settings.MaxIter; iter++) {

    // Compute a few quantities
    normGrad = vector_infnorm(&solver->grad);
    df = fabs(fOld-f);

    // Save old values
    fOld = f;
    matrix_copy(x, &solver->xOld);
    matrix_copy(&solver->grad, &solver->gradOld);
    matrix_copy(&solver->hessian, &solver->hessianOld);

    // Test for convergence
    if (normGrad < solver->settings.TolFun && positiveDefinite)
      break;
    if (normDx < 0.9f*delta && ratio > 0.25f && df < solver->settings.TolFun*(1.0f+fabs(fOld)))
      break;
    if (normDx < solver->settings.TolX)
      break;
    
    // Get new x
    value = step(solver, x, delta, computeZ, &positiveDefinite);
    normDx = vector_norm(&solver->dx);
    matrix_add(x, &solver->dx, x);

    // Evaluate
    f = solver->functions.f(x, &solver->grad, &solver->hessian, solver->functions.p);
    
    // Update
    ratio = (f - fOld)/value;
    if (ratio >= 0.75f && normDx >= 0.9f*delta)
      delta = MIN(delbnd, 2.0f*delta);
    else if (ratio <= 0.25f)
      delta = MIN(normDx/4.0f, delta/4.0f);
        
    // Is the new point better?
    if (fOld <= f) {
      // No, revert
      f = fOld;
      matrix_copy(&solver->xOld, x);
      matrix_copy(&solver->gradOld, &solver->grad);
      matrix_copy(&solver->hessianOld, &solver->hessian);
      computeZ = 0;
    } else {
      // Yes, we can reset Z
      computeZ = 1;
    }

    // Print
    if (solver->settings.verbose)
      printf("#%d           %10f\n", iter, f);
  }
}

/****************************************************************
 * The following function computes the CG-step                  *
 ****************************************************************/

static float step(pcgtr_solver_t *solver, matrix_t *x, float delta, int computeZ, int *positiveDefinite) {
  float normV1, normV2, V1dotV2;
  float value;

  *positiveDefinite = 1;

  // Preconditioning (finding Z in 2-dimensions - v1, v2)
  if (computeZ) {
    matrix_zero(&solver->Z);

    // Preconditioning loop
    precondition(solver, positiveDefinite);

    // Set Z
    normV1 = vector_norm(&solver->v1);
    if (normV1 > solver->settings.eps) {
      matrix_setcol_real(&solver->v1, 1.0f/normV1, 0, &solver->Z);
    }

    if (!*positiveDefinite) {
      vector_sign(&solver->grad, &solver->v2);
    } else {
      matrix_copy(&solver->grad, &solver->v2);
    }
    normV2 = vector_norm(&solver->v2);
    if (normV2 > solver->settings.eps) {
      vector_mul_real(&solver->v2, 1.0f/normV2, &solver->v2);
    }

    V1dotV2 = vector_scalar_product(&solver->v1, &solver->v2);
    vector_sub_real_vector(&solver->v2, V1dotV2, &solver->v1, &solver->v2);

    normV2 = vector_norm(&solver->v2);
    if (normV2 > solver->settings.eps) {
      matrix_setcol_real(&solver->v2, 1.0f/normV2, 1, &solver->Z);
    }
  }

  // Reduce to the chosen subspace
  transpose_matrix_mul_vector(&solver->Z, &solver->grad, &solver->trustG);
  transpose_matrix_mul(&solver->Z, solver->functions.HessianMultiplication(&solver->hessian, &solver->Z, &solver->A, solver->functions.p), &solver->trustH);
  // Determine 2D trust region solution: min{g'*s + 0.5*s'*H*s} such that ||s|| <= delta
  trustregion(solver, delta);
  // Transform back to nominal space
  matrix_mul(&solver->Z, &solver->trustS, &solver->dx);
  // Value of the TR solution
  value = vector_scalar_product(&solver->trustG, &solver->trustS) + 0.5f*vector_scalar_product(&solver->trustS, matrix_mul(&solver->trustH, &solver->trustS, &solver->B));

  return value;
}

/**************************************************************
 * The following function implements the preconditioning loop *
 * for the conjugate gradient                                 *
 **************************************************************/

static void precondition(pcgtr_solver_t *solver, int *positiveDefinite) {
  int k;
  float TolStop, pTAp, numerator, denominator, beta, alpha, normP;
  matrix_t *x = &solver->v1;
  matrix_t *r = &solver->rpcg;
  matrix_t *z = &solver->zpcg;
  matrix_t *p = &solver->ppcg;
  matrix_t *Ap = &solver->Appcg;
  matrix_t *M = &solver->Mpcg;

  /* We approximate the hessian with a diagonal matrix such that H ~ R'*R.
     M is the inverse of the approximated H (inv(H) ~ inv(R)*inv(R').
     Note that M here is the not the M matrix but its diagonal. */
  solver->functions.DiagonalInverseHessianApproximation(&solver->hessian, M, solver->functions.p);

  /* Preconditioned conjugate gradient method using Fletcher-Reeves
     (http://en.wikipedia.org/wiki/Conjugate_gradient_method) */
  
  // Setup
  *positiveDefinite = 1;
  vector_zero(x);
  vector_neg(&solver->grad, r);
  vector_mul_elements(M, r, z);
  TolStop = 0.1f*vector_norm(z);

  // Helper (beta is z'*r at k+1 divided by z'r at k, beta = numerator/denominator)
  numerator = vector_scalar_product(r, z);
  denominator = 0.0f;

  for (k = 0; k < solver->settings.MaxPCGIter; k++) {
    if (k == 0) {
      vector_copy(z, p);
    } else {
      beta = numerator/denominator;
      vector_add_real_vector(z, beta, p, p);
    }

    solver->functions.HessianMultiplication(&solver->hessian, p, Ap, solver->functions.p);
    pTAp = vector_scalar_product(p, Ap);

    if (pTAp <= 0.0f) {
      normP = vector_norm(p);
      if (normP < solver->settings.eps) {
        vector_copy(p, x);
      } else {
        vector_mul_real(p, 1.0f/normP, x);
      }
      *positiveDefinite = 0;
      break;
    } else {
      alpha = numerator/pTAp;
      vector_add_real_vector(x, alpha, p, x);
      vector_sub_real_vector(r, alpha, Ap, r);
    }

    vector_mul_elements(M, r, z);
    if (vector_norm(p) < TolStop) {
      break;
    }
    denominator = numerator;
    numerator = vector_scalar_product(r, z);
  }
}

/*************************************************************
 * The following function implements the trust region search *
 *************************************************************/

static void trustregion(pcgtr_solver_t *solver, float delta) {
  // We avoid using vectors/matrix operation as they are slower when n is small (2 in this case)
  float alpha[2];
  float coefficients[2];
  float *eigenvalues;
  float *eigenvectors;
  float smallestEigenvalue;
  float lambda, lambdaInit;
  float normS, value;
  int smallestIndex;
  int sig;
  float *s = solver->trustS.v;
  float *g = solver->trustG.v;
  matrix_t *H = &solver->trustH;

  // Get eigenvalues/eigenvectors of H
  eigenvalues = solver->eigenvaluesTR.v;
  eigenvectors = solver->eigenvectorsTR.v;
  matrix_eig(H, &solver->eigenvectorsTR, &solver->eigenvaluesTR);

  // Get smallest eigenvalue index
  smallestIndex = 0;
  smallestEigenvalue = eigenvalues[0];
  if (eigenvalues[1] < smallestEigenvalue) {
    smallestIndex = 1;
    smallestEigenvalue = eigenvalues[1];
  }
  
  // Set alpha = -eigenvectors'*g
  alpha[0] = -eigenvectors[0]*g[0] - eigenvectors[2]*g[1];
  alpha[1] = -eigenvectors[1]*g[0] - eigenvectors[3]*g[1];

  // Sign of alpha
  sig = (alpha[smallestIndex] >= 0.0f) - (alpha[smallestIndex] < 0.0f);

  // initialize coefficients
  coefficients[0] = 0.0f;
  coefficients[1] = 0.0f;

  if (smallestEigenvalue > 0.0f) {
    // s = eigenvectors*alpha ./ eigenvalues
    coefficients[0] = alpha[0] / eigenvalues[0];
    coefficients[1] = alpha[1] / eigenvalues[1];
    s[0] = eigenvectors[0]*coefficients[0] + eigenvectors[1]*coefficients[1];
    s[1] = eigenvectors[2]*coefficients[0] + eigenvectors[3]*coefficients[1];
    normS = sqrtf(s[0]*s[0] + s[1]*s[1]);
    if (normS <= 1.2f*delta) {
      // Good value for s
      return;
    } else {
      // Bad value for s
      lambdaInit = 0.0f;
    }
  } else {
    // H is negative definite
    lambdaInit = -smallestEigenvalue;
  }

  if (secular_equation(solver, lambdaInit, eigenvalues, alpha, delta) > 0.0f) {
    lambda = rfzero(solver, lambdaInit, eigenvalues, alpha, delta, &value);
    if (value < solver->settings.TolZero) {
      if (get_solution(solver, lambda, eigenvalues, eigenvectors, alpha, s, 0)) {
        lambda = -smallestEigenvalue;
      } else {
        normS = sqrtf(s[0]*s[0] + s[1]*s[1]);
        if (normS > 1.2f*delta || normS < 0.8f*delta) {
          lambda = -smallestEigenvalue;
        } else {
          // We are happy
          return;
        }
      }
    } else {
      lambda = -smallestEigenvalue;
    }
  } else {
    lambda = -smallestEigenvalue;
  }
  
  // We had bad values until now (limit the damage -- use failsafe mode)
  get_solution(solver, lambda, eigenvalues, eigenvectors, alpha, s, 1);
  normS = sqrtf(s[0]*s[0] + s[1]*s[1]);

  if (normS < 0.8f*delta) {
    // Norm is too small, artifically increase it
    float beta = sqrtf(delta*delta + normS*normS);
    s[0] += beta*sig*eigenvectors[smallestIndex];
    s[1] += beta*sig*eigenvectors[smallestIndex+2];
  } else if (normS > 1.2f*delta) {
    // Norm is too big, correct alpha and find a new lambda
    float w[2];
    w[0] = eigenvalues[0] + lambda;
    w[1] = eigenvalues[1] + lambda;
    if (fabs(w[0]) < solver->settings.eps) alpha[0] = 0.0f;
    if (fabs(w[1]) < solver->settings.eps) alpha[1] = 0.0f;
    lambda = rfzero(solver, lambdaInit, eigenvalues, alpha, delta, &value);
    // Get the new solution (failsafe)
    get_solution(solver, lambda, eigenvalues, eigenvectors, alpha, s, 1);
  }
}

static int get_solution(pcgtr_solver_t *solver, float lambda, float *eigenvalues, float *eigenvectors, float *alpha, float *s, int failsafe) {
  float w[2];
  float coefficients[2];

  // Compute the solution and verify that it is decent
  w[0] = eigenvalues[0] + lambda;
  w[1] = eigenvalues[1] + lambda;

  // Avoid division by zero (alpha[0]/w[0])
  if (fabs(w[0]) < solver->settings.eps) {
    if (fabs(alpha[0]) < solver->settings.eps || failsafe) coefficients[0] = 0.0f;
    else return 1;
  } else coefficients[0] = alpha[0]/w[0];

  // Avoid division by zero (alpha[1]/w[1])
  if (fabs(w[1]) < solver->settings.eps) {
    if (fabs(alpha[1]) < solver->settings.eps || failsafe) coefficients[1] = 0.0f;
    else return 1;
  } else coefficients[1] = alpha[1]/w[1];

  // Compute s
  s[0] = eigenvectors[0]*coefficients[0] + eigenvectors[1]*coefficients[1];
  s[1] = eigenvectors[2]*coefficients[0] + eigenvectors[3]*coefficients[1];
  return 0;
}

static float secular_equation(pcgtr_solver_t *solver, float lambda, float *eigenvalues, float *alpha, float delta) {
  float nlambda[2];

  nlambda[0] = eigenvalues[0] + lambda;
  nlambda[1] = eigenvalues[1] + lambda;
  if (fabs(nlambda[0]) < solver->settings.eps || fabs(nlambda[1]) < solver->settings.eps) return 1.0f/delta;

  nlambda[0] = alpha[0] / nlambda[0];
  nlambda[1] = alpha[1] / nlambda[1];
  return 1.0f/delta - 1.0f/sqrtf(nlambda[0]*nlambda[0] + nlambda[1]*nlambda[1]);
}

static float rfzero(pcgtr_solver_t *solver, float x, float *eigenvalues, float *alpha, float delta, float *value) {
  int iter = 0;
  float dx, a, b, c, fa, fb, fc, m, s, p, q, r, tol;
  float d = 0.0f;
  float e = 0.0f;

  // Finds the zero to the right of x (of the secular_equation)
    
  // Initialization
  dx = fabs(x)/2.0f;
  if (x == 0.0f) dx = 0.5f;
    
  a = x;
  c = x;
  fa = secular_equation(solver, a, eigenvalues, alpha, delta);
  iter++;
        
  b = x + dx;
  fb = secular_equation(solver, b, eigenvalues, alpha, delta);
  iter++;

  // Find the change of sign
  while ((fa > 0.0f) == (fb > 0.0f)) {
    // Increase dx
    dx *= 2.0f;
    b = x + dx;
    fb = secular_equation(solver, b, eigenvalues, alpha, delta);
    iter++;
    if (iter > solver->settings.MaxZeroIter) break;
  }
        
  fc = fb;
  while (1) {
    // Ensure that b is the best result (a is the previous value of
    // b and c is on the opposite side of the zero from b).
    if ((fb > 0.0f) == (fc > 0.0f)) {
      c = a;
      fc = fa;
      d = b-a;
      e = d;
    }
    if (fabs(fc) < fabs(fb)) {
      // Swap b and c
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    if (iter >= solver->settings.MaxZeroIter) break;
    
    m = 0.5f*(c-b);
    tol = 2.0f*solver->settings.TolZero*MAX(fabs(b), 1.0f);
    if (fabs(m) <= tol || fb == 0.0f)
      break;
            
    // Bisect
    if (fabs(e) < tol || fabs(fa) <= fabs(fb)) {
      d = m;
      e = m;
    } else {
      s = fb/fa;
      if (a == c) {
        // Linear interpolation
        p = 2.0f*m*s;
        q = 1.0f - s;
      } else {
        // Quadratic interpolation
        q = fa/fc;
        r = fb/fc;
        p = s*(2.0f*m*q*(q - r) - (b - a)*(r - 1.0f));
        q = (q - 1.0f)*(r - 1.0f)*(s - 1.0f);
      }
                
      if (p > 0) q = -q;
      else p = -p;
                
      // Is interpolated point acceptable?
      if ((2.0f*p < 3.0f*m*q - fabs(tol*q)) && (p < fabs(0.5f*e*q))) {
        e = d;
        d = p/q;
      } else {
        d = m; 
        e = m;
      }
    }
            
    // Next point
    a = b;
    fa = fb;
    if (fabs(d) > tol)
      b = b + d;
    else if (b > c)
      b = b - tol;
    else
      b = b + tol;
    fb = secular_equation(solver, b, eigenvalues, alpha, delta);
    iter++;
  }

  *value = fb;
  return b;
}

/***********************************************************
 * The following function approximates the hessian inverse *
 ***********************************************************/

static matrix_t *DefaultDiagonalInverseHessianApproximation(matrix_t *a, matrix_t *z, void *p) {
  int i, j;
  int m = a->m;
  int n = a->n;
  float s, v;

  for (i = 0; i < m; i++) {
    s = 0.0f;
    for (j = 0; j < n; j++) {
      v = a->v[i*n+j];
      s += v*v;
    }
    z->v[i] = 1.0f/MAX(sqrtf(s), 0.01f);
  }

  return z;
}

/******************************************************************
 * The following function computes a normal matrix multiplication *
 ******************************************************************/

static matrix_t *DefaultHessianMultiplication(matrix_t *a, matrix_t *b, matrix_t *z, void *p) {
  return matrix_mul(a, b, z);
}
