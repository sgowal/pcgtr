/**
 * \brief A stripped matrix library.
 * This file provides a simple matrix library.
 * All the function provided do NOT make any
 * verification on either the size of the matrices
 * or whether auxilary functions (i.e. malloc) succeed.
 */

#ifndef MATRIX_H
#define MATRIX_H

/**
 * \brief A simple matrix structure.
 */
typedef struct {
  int m;
  int n;
  float *v;
} matrix_t;

/**
 * \brief Creates a matrix with the number of rows and columns indicated.
 * @param m a pointer to the matrix to create
 * @param rows the number of rows
 * @param cols the number of columns
 * @return m (the same pointer given as input)
 */
matrix_t *matrix_init(matrix_t *m, int rows, int cols);

/**
 * \brief Destroy (deallocates the matrix m).
 * @param m a pointer to the matrix to destroy
 */
void matrix_destroy(matrix_t *m);

/**
 * \brief Returns the index in the value array that corresponds to row and col.
 * @param row 0-index of the row
 * @param col 0-index of the vol
 * @returns index in the value array that corresponds to row and col.
 */
int matrix_index(matrix_t *m, int row, int col);

/**
 * \brief Sets all elements of the matrix to 0
 * @param a a pointer to the matrix
 * @return a (the same pointer given as input)
 */
matrix_t *matrix_zero(matrix_t *a);

/**
 * \brief Sets all elements of the vector to 0
 * @param a a pointer to the vector
 * @return a (the same pointer given as input)
 */
matrix_t *vector_zero(matrix_t *a);

/**
 * \brief Copies the matrix a into z.
 * @param a a pointer to the matrix a
 * @param z a pointer to the matrix z
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_copy(matrix_t *a, matrix_t *z);

/**
 * \brief Copies the vector a into z.
 * @param a a pointer to the vector a
 * @param z a pointer to the vector z
 * @return z (the same pointer given as input)
 */
matrix_t *vector_copy(matrix_t *a, matrix_t *z);

/**
 * \brief Copies the vector x into the ith row of matrix z.
 * @param x a pointer to the vector x
 * @param i the row index
 * @param z a pointer to the matrix z
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_setrow(matrix_t *x, int i, matrix_t *z);

/**
 * \brief Copies the vector x into the jth column of matrix z.
 * @param x a pointer to the vector x
 * @param j the column index
 * @param z a pointer to the matrix z
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_setcol(matrix_t *x, int j, matrix_t *z);

/**
 * \brief Copies the vector x*a into the jth column of matrix z.
 * @param x a pointer to the vector x
 * @param j the column index
 * @param a the real number
 * @param z a pointer to the matrix z
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_setcol_real(matrix_t *x, float a, int j, matrix_t *z);

/**
 * \brief Computes the 2-norm of the vector x
 * @param x a pointer to the vector x
 * @return the 2-norm of x
 */
float vector_norm(matrix_t *x);

/**
 * \brief Computes the squared 2-norm of the vector x
 * @param x a pointer to the vector x
 * @return the squared 2-norm of x
 */
float vector_2norm(matrix_t *x);

/**
 * \brief Computes the norm of the vector x
 * @param x a pointer to the vector x
 * @return the infinity-norm of x
 */
float vector_infnorm(matrix_t *x);

/**
 * \brief Computes the sign of each element of the vector x and store them in z.
 * The sign of an element is -1 if negative, 1 if positive and 0 if equal to 0.
 * @param x a pointer to the vector x
 * @return z (the same pointer given as input)
 */
matrix_t *vector_sign(matrix_t *x, matrix_t *z);

/**
 * \brief Multiplies a vector with a real number and store in z.
 * @param x a pointer to the vector x
 * @param a the real number
 * @param z a pointer to the matrix z (this can be x)
 * @return z (the same pointer given as input)
 */
matrix_t *vector_mul_real(matrix_t *x, float a, matrix_t *z);

/**
 * \brief Does the scalar product of x and y.
 * @param x a pointer to the vector x
 * @param y a pointer to the vector y
 * @return the scalar product
 */
float vector_scalar_product(matrix_t *x, matrix_t *y);

/**
 * \brief Adds two vectors a and c*b (this 1.5x faster than the matrix_add(matrix_mul_real) version).
 * @param a a pointer to vector a
 * @param b a pointer to vector b
 * @param c the real number
 * @param z a pointer to the resulting vector (this can be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *vector_add_real_vector(matrix_t *a, float c, matrix_t *b, matrix_t *z);

/**
 * \brief Substracts two vectors a and c*b (this 1.5x faster than the matrix_add(matrix_mul_real) version).
 * @param a a pointer to vector a
 * @param b a pointer to vector b
 * @param c the real number
 * @param z a pointer to the resulting vector (this can be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *vector_sub_real_vector(matrix_t *a, float c, matrix_t *b, matrix_t *z);

/**
 * \brief Multiplies element-wise a vector with a with b.
 * @param a a pointer to the vector a
 * @param b a pointer to the vector b
 * @param z a pointer to the vector z (this can be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *vector_mul_elements(matrix_t *a, matrix_t *b, matrix_t *z);

/**
 * \brief Multiplies a by -1.
 * @param a a pointer to the vector a
 * @param z a pointer to the vector z (this can be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *vector_neg(matrix_t *a, matrix_t *z);

/**
 * \brief Multiplies two matrices a and b.
 * @param a a pointer to matrix a
 * @param b a pointer to matrix b
 * @param z a pointer to the resulting matrix (this cannot be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_mul(matrix_t *a, matrix_t *b, matrix_t *z);

/**
 * \brief Multiplies a matrix a and a vector b.
 * @param a a pointer to matrix a
 * @param b a pointer to vector b
 * @param z a pointer to the resulting vector (this cannot be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_mul_vector(matrix_t *a, matrix_t *b, matrix_t *z);

/**
 * \brief Multiplies two matrices a and b^T.
 * @param a a pointer to matrix a
 * @param b a pointer to matrix b
 * @param z a pointer to the resulting matrix (this cannot be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_mul_transpose(matrix_t *a, matrix_t *b, matrix_t *z);

/**
 * \brief Multiplies two matrices a^T and b.
 * @param a a pointer to matrix a
 * @param b a pointer to matrix b
 * @param z a pointer to the resulting matrix (this cannot be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *transpose_matrix_mul(matrix_t *a, matrix_t *b, matrix_t *z);

/**
 * \brief Multiplies a matrix a^T and a vector b.
 * @param a a pointer to matrix a
 * @param b a pointer to vector b
 * @param z a pointer to the resulting vector (this cannot be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *transpose_matrix_mul_vector(matrix_t *a, matrix_t *b, matrix_t *z);

/**
 * \brief Multiplies a matrice with real number.
 * @param a a pointer to matrix a
 * @param b the real number
 * @param z a pointer to the resulting matrix (this can be a)
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_mul_real(matrix_t *a, float b, matrix_t *z);

/**
 * \brief Multiplies two matrices a and c*b (this 1.5x faster than the matrix_mul(matrix_mul_real) version).
 * @param a a pointer to matrix a
 * @param b a pointer to matrix b
 * @param c the real number
 * @param z a pointer to the resulting matrix (this cannot be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_mul_real_matrix(matrix_t *a, float c, matrix_t *b, matrix_t *z);

/**
 * \brief Multiplies two matrices a and c*b^T.
 * @param a a pointer to matrix a
 * @param b a pointer to matrix b
 * @param z a pointer to the resulting matrix (this cannot be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_mul_real_transpose(matrix_t *a, float c, matrix_t *b, matrix_t *z);

/**
 * \brief Multiplies two matrices a^T and c*b.
 * @param a a pointer to matrix a
 * @param b a pointer to matrix b
 * @param z a pointer to the resulting matrix (this cannot be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *transpose_matrix_mul_real_matrix(matrix_t *a, float c, matrix_t *b, matrix_t *z);

/**
 * \brief Multiplies two matrices a and -b (this 1.5x faster than the matrix_mul(matrix_mul_real) version).
 * @param a a pointer to matrix a
 * @param b a pointer to matrix b
 * @param z a pointer to the resulting matrix (this cannot be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_mul_neg_matrix(matrix_t *a, matrix_t *b, matrix_t *z);

/**
 * \brief Adds two matrices a and b.
 * @param a a pointer to matrix a
 * @param b a pointer to matrix b
 * @param z a pointer to the resulting matrix (this can be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_add(matrix_t *a, matrix_t *b, matrix_t *z);

/**
 * \brief Adds two matrices a and c*b (this 1.5x faster than the matrix_add(matrix_mul_real) version).
 * @param a a pointer to matrix a
 * @param b a pointer to matrix b
 * @param c the real number
 * @param z a pointer to the resulting matrix (this can be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_add_real_matrix(matrix_t *a, float c, matrix_t *b, matrix_t *z);

/**
 * \brief Substract two matrices a and b.
 * @param a a pointer to matrix a
 * @param b a pointer to matrix b
 * @param z a pointer to the resulting matrix (this can be a or b)
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_sub(matrix_t *a, matrix_t *b, matrix_t *z);

/**
 * \brief Squares all elements of a matrix a and sums all rows into z.
 * @param a a pointer to matrix a
 * @param z a pointer to the resulting vector
 * @return z (the same pointer given as input)
 */
matrix_t *matrix_square_sum(matrix_t *a, matrix_t *z);

/**
 * \brief Computes the eigenvectors and eigenvalues of a 2x2 matrix a
 * @param a a pointer to matrix a (2x2)
 * @param V a pointer to the resulting eigenvector matrix (2x2)
 * @param D a pointer to the resulting eigenvalues vector (2x1)
 */
void matrix_eig(matrix_t *a, matrix_t *V, matrix_t *D);

/**
 * \brief Prints a matrix with the name specified.
 * @param a a pointer to the matrix
 * @param name the name of the matrix a
 */
void matrix_print(matrix_t *a, char *name);
void array_print(float *a, char *name, int m, int n);

/**
 * \brief Prints a matrix transposed with the name specified.
 * @param a a pointer to the matrix
 * @param name the name of the matrix a
 */
void transpose_matrix_print(matrix_t *a, char *name);

#endif
