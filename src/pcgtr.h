#ifndef PCGTR_H
#define PCGTR_H

#include "matrix.h"

typedef struct {
  /** Number of iterations */
  int MaxIter;
  /** Number of PCG iterations (each time the Newton direction needs to be found) */
  int MaxPCGIter;
  /** Number of iterations for the trust-region solution */
  int MaxZeroIter;
  /** Tolerance on the function for declaring convergence */
  float TolFun;
  /** Tolerance on the x-values */
  float TolX;
  /** Tolerance on the trust-region solution */
  float TolZero;
  /** Value of an epsilon value (to avoid division by zero) */
  float eps;
  /** Verbosity */
  int verbose;
} pcgtr_settings_t;

typedef struct {
  /** Function to optimize. It compute the Jacobian and the Hessian as well */
  float (*f)(matrix_t *x, matrix_t *J, matrix_t *H, void *p);
  /** Function that multiplies the Hessian H with another matrix Y and returns the result in Z.
      It is initially set to matrix_mul. */
  matrix_t* (*HessianMultiplication)(matrix_t *H, matrix_t *Y, matrix_t *Z, void *p);
  /** Function that set the vector M to the diagonal of the approximated Hessian H.
      It is set to an internal function that returns a vector with its i-th element
      equal to the inverse squared sum of the element of the Hessian on the i-th row. */
  matrix_t* (*DiagonalInverseHessianApproximation)(matrix_t *H, matrix_t *M, void *p);
  /** Optional data to give as input to those functions */
  void *p;
} pcgtr_functions_t;

typedef struct {
  /** Settings */
  pcgtr_settings_t settings;

  /** Functions */
  pcgtr_functions_t functions;

  /** Storage -- DO NOT USE */
  matrix_t v1;
  matrix_t v2;
  matrix_t trustG;
  matrix_t trustH;
  matrix_t A;
  matrix_t B;
  matrix_t trustS;
  matrix_t rpcg;
  matrix_t zpcg;
  matrix_t ppcg;
  matrix_t Appcg;
  matrix_t Mpcg;
  matrix_t eigenvaluesTR;
  matrix_t eigenvectorsTR;
  matrix_t grad;
  matrix_t hessian;
  matrix_t Z;
  matrix_t xOld;
  matrix_t gradOld;
  matrix_t hessianOld;
  matrix_t dx;
} pcgtr_solver_t;

/**
 * \brief Initializes the cgsolver.
 * This function allocates the necessary memory.
 * @param n the number of variables that will be optimized
 * @param sparse if true will not initialize the Hessian matrix.
 * @return 0 on success
 */
int pcgtr_solver_init(pcgtr_solver_t *solver, int n, int sparse);

/**
 * \brief Stops the cgsolver.
 * This function deallocates the memory.
 */
void pcgtr_solver_destroy(pcgtr_solver_t *solver);

/**
 * \brief Replaces usual NxN Hessian matrix by a row x col matrix.
 * This function allocates the Hessian matrix as a row x col matrix instead.
 * It should only use in sparse mode when HessianMultiplication and
 * DiagonalInverseHessianApproximation functions are manually defined.
 * @param rows the number of rows
 * @param cols the number of columns
 */
void pcgtr_solver_halloc(pcgtr_solver_t *solver, int row, int col);

/**
 * \brief Solve the unconstrained minimization problem min{f(x)}
 * This function uses the Conjugate Gradient method to minimize a cost
 * function f. The function f is responsible to generate both the gradient
 * of the cost J and the hessian H.
 * @param f a cost function
 * @param x both the initial solution and final solution
 * @param hess_mult a function that multiplies the Hessian with another matrix Y (H*Y = Z), if NULL is passed the standard matrix multiplication is used.
 * @param hess_approx a function that provides the inverse square root of the sum of the squares of each row of the hessian (1./sum(H.^2,2) = M), if NULL is passed the standard formula is used.
 * @param p a pointer to some user data needed by f (as fourth argument)
 */
void pcgtr_solver_solve(pcgtr_solver_t *solver, matrix_t *x /**< [in,out] */);

#endif
