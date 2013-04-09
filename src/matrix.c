#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"

matrix_t *matrix_init(matrix_t *m, int rows, int cols) {
  m->m = rows;
  m->n = cols;
  m->v = (float *)calloc(sizeof(float), rows*cols);
  if (!m->v)
    return NULL;
  return m;
}

void matrix_destroy(matrix_t *m) {
  free(m->v);
}

int matrix_index(matrix_t *m, int row, int col) {
  return row*m->n + col;
}

matrix_t *matrix_zero(matrix_t *a) {
  int s = a->m*a->n;
  memset(a->v, 0, s*sizeof(float));
  return a;
}

matrix_t *vector_zero(matrix_t *a) {
  int s = a->m;
  memset(a->v, 0, s*sizeof(float));
  return a;
}

matrix_t *matrix_copy(matrix_t *a, matrix_t *z) {
  int s = a->m*a->n;
  memcpy(z->v, a->v, s*sizeof(float));
  return z;
}

matrix_t *vector_copy(matrix_t *a, matrix_t *z) {
  int s = a->m;
  memcpy(z->v, a->v, s*sizeof(float));
  return z;
}

matrix_t *matrix_setrow(matrix_t *x, int i, matrix_t *z) {
  int n = z->n;
  memcpy(z->v + i*n, x->v, n*sizeof(float));
  return z;
}

matrix_t *matrix_setcol(matrix_t *x, int j, matrix_t *z) {
  int i;
  int n = z->n;
  int m = z->m;

  for (i = 0; i < m; i++) {
    z->v[j] = x->v[i];
    j += n;
  }

  return z;
}

matrix_t *matrix_setcol_real(matrix_t *x, float a, int j, matrix_t *z) {
  int i;
  int n = z->n;
  int m = z->m;

  for (i = 0; i < m; i++) {
    z->v[j] = a*x->v[i];
    j += n;
  }

  return z;
}

float vector_norm(matrix_t *x) {
  int i;
  int n = x->m;
  float norm = 0.0f;
  
  for (i = 0; i < n; i++) {
    norm += x->v[i]*x->v[i];
  }

  return sqrtf(norm);
}

float vector_2norm(matrix_t *x) {
  int i;
  int n = x->m;
  float norm = 0.0f;
  
  for (i = 0; i < n; i++) {
    norm += x->v[i]*x->v[i];
  }

  return norm;
}

float vector_infnorm(matrix_t *x) {
  int i;
  int n = x->m;
  float norm = -1.0f;
  float v;
  
  for (i = 0; i < n; i++) {
    v = fabs(x->v[i]);
    if (v > norm) norm = v;
  }

  return norm;
}

matrix_t *vector_sign(matrix_t *x, matrix_t *z) {
  int i;
  int n = x->m;

  for (i = 0; i < n; i++) {
    z->v[i] = (float)((x->v[i] > 0.0f) - (x->v[i] < 0.0f));
  }

  return z;
}

float vector_scalar_product(matrix_t *x, matrix_t *y) {
  int i;
  int n = x->m;
  float s = 0.0f;

  for (i = 0; i < n; i++) {
    s += x->v[i]*y->v[i];
  }

  return s;
}

matrix_t *matrix_mul(matrix_t *a, matrix_t *b, matrix_t *z) {  
  int i, j, k;
  int m = a->m;
  int n = a->n;
  int o = b->n;
  float s;

  for (i = 0; i < m; i++) {
    for (j = 0; j < o; j++) {
      s = 0.0f;
      for (k = 0; k < n; k++) {
        s += a->v[i*n+k]*b->v[k*o+j];
      }
      z->v[i*o+j] = s;
    }
  }

  return z;
}

matrix_t *matrix_mul_vector(matrix_t *a, matrix_t *b, matrix_t *z) {  
  int i, k;
  int m = a->m;
  int n = a->n;
  float s;

  for (i = 0; i < m; i++) {
    s = 0.0f;
    for (k = 0; k < n; k++) {
      s += a->v[i*n+k]*b->v[k];
    }
    z->v[i] = s;
  }

  return z;
}

matrix_t *matrix_mul_transpose(matrix_t *a, matrix_t *b, matrix_t *z) {
  int i, j, k;
  int m = a->m;
  int n = a->n;
  int o = b->m;
  float s;

  for (i = 0; i < m; i++) {
    for (j = 0; j < o; j++) {
      s = 0.0f;
      for (k = 0; k < n; k++) {
        s += a->v[i*n+k]*b->v[j*o+k];
      }
      z->v[i*o+j] = s;
    }
  }

  return z;
}

matrix_t *transpose_matrix_mul(matrix_t *a, matrix_t *b, matrix_t *z) {
  int i, j, k;
  int m = a->m;
  int n = a->n;
  int o = b->n;
  float s;

  for (i = 0; i < n; i++) {
    for (j = 0; j < o; j++) {
      s = 0.0f;
      for (k = 0; k < m; k++) {
        s += a->v[k*n+i]*b->v[k*o+j];
      }
      z->v[i*o+j] = s;
    }
  }

  return z;
}

matrix_t *transpose_matrix_mul_vector(matrix_t *a, matrix_t *b, matrix_t *z) {
  int i, k;
  int m = a->m;
  int n = a->n;
  float s;

  for (i = 0; i < n; i++) {
    s = 0.0f;
    for (k = 0; k < m; k++) {
      s += a->v[k*n+i]*b->v[k];
    }
    z->v[i] = s;
  }

  return z;
}

matrix_t *matrix_mul_real(matrix_t *a, float b, matrix_t *z) {
  int i;
  int n = a->m*a->n;

  for (i = 0; i < n; i++) {
    z->v[i] = b*a->v[i];
  }

  return z;
}

matrix_t *vector_mul_real(matrix_t *x, float a, matrix_t *z) {
  int i;
  int n = x->m;

  for (i = 0; i < n; i++) {
    z->v[i] = a*x->v[i];
  }

  return z;
}

matrix_t *vector_add_real_vector(matrix_t *a, float c, matrix_t *b, matrix_t *z) {
  int i;
  int n = a->m;

  for (i = 0; i < n; i++)
    z->v[i] = a->v[i]+c*b->v[i];

  return z;
}

matrix_t *vector_sub_real_vector(matrix_t *a, float c, matrix_t *b, matrix_t *z) {
  int i;
  int n = a->m;

  for (i = 0; i < n; i++)
    z->v[i] = a->v[i]-c*b->v[i];

  return z;
}

matrix_t *vector_mul_elements(matrix_t *a, matrix_t *b, matrix_t *z) {
  int i;
  int n = a->m;

  for (i = 0; i < n; i++)
    z->v[i] = a->v[i]*b->v[i];

  return z;
}

matrix_t *vector_neg(matrix_t *a, matrix_t *z) {
  int i;
  int n = a->m;

  for (i = 0; i < n; i++)
    z->v[i] = -a->v[i];

  return z;
}

matrix_t *matrix_mul_real_matrix(matrix_t *a, float c, matrix_t *b, matrix_t *z) {
  int i, j, k;
  int m = a->m;
  int n = a->n;
  int o = b->n;
  float s;

  for (i = 0; i < m; i++) {
    for (j = 0; j < o; j++) {
      s = 0.0f;
      for (k = 0; k < n; k++) {
        s += c*a->v[i*n+k]*b->v[k*o+j];
      }
      z->v[i*o+j] = s;
    }
  }

  return z;
}

matrix_t *matrix_mul_real_transpose(matrix_t *a, float c, matrix_t *b, matrix_t *z) {
  int i, j, k;
  int m = a->m;
  int n = a->n;
  int o = b->m;
  float s;

  for (i = 0; i < m; i++) {
    for (j = 0; j < o; j++) {
      s = 0.0f;
      for (k = 0; k < n; k++) {
        s += c*a->v[i*n+k]*b->v[j*o+k];
      }
      z->v[i*o+j] = s;
    }
  }

  return z;
}

matrix_t *transpose_matrix_mul_real_matrix(matrix_t *a, float c, matrix_t *b, matrix_t *z) {
  int i, j, k;
  int m = a->m;
  int n = a->n;
  int o = b->n;
  float s;

  for (i = 0; i < n; i++) {
    for (j = 0; j < o; j++) {
      s = 0.0f;
      for (k = 0; k < m; k++) {
        s += c*a->v[k*n+i]*b->v[k*o+j];
      }
      z->v[i*o+j] = s;
    }
  }

  return z;
}

matrix_t *matrix_mul_neg_matrix(matrix_t *a, matrix_t *b, matrix_t *z) {
  int i, j, k;
  int m = a->m;
  int n = a->n;
  int o = b->n;
  float s;

  for (i = 0; i < m; i++) {
    for (j = 0; j < o; j++) {
      s = 0.0f;
      for (k = 0; k < n; k++) {
        s += -a->v[i*n+k]*b->v[k*o+j];
      }
      z->v[i*o+j] = s;
    }
  }

  return z;
}

matrix_t *matrix_add(matrix_t *a, matrix_t *b, matrix_t *z) {
  int i;
  int n = a->m*a->n;

  for (i = 0; i < n; i++)
    z->v[i] = a->v[i]+b->v[i];

  return z;
}

matrix_t *matrix_add_real_matrix(matrix_t *a, float c, matrix_t *b, matrix_t *z) {
  int i;
  int n = a->m*a->n;

  for (i = 0; i < n; i++)
    z->v[i] = a->v[i]+c*b->v[i];

  return z;
}

matrix_t *matrix_sub(matrix_t *a, matrix_t *b, matrix_t *z) {
  int i;
  int n = a->m*a->n;

  for (i = 0; i < n; i++)
    z->v[i] = a->v[i]-b->v[i];

  return z;
}

matrix_t *matrix_square_sum(matrix_t *a, matrix_t *z) {
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
    z->v[i] = s;
  }

  return z;
}

void matrix_eig(matrix_t *A, matrix_t *V, matrix_t *D) {
  float a = A->v[0];
  float b = A->v[1];
  float c = A->v[2];
  float d = A->v[3];
  
  float v1x, v1y, v2x, v2y, l1, l2;
  if (c < 0.00001f && c > -0.00001f) {
    if (a - d < 0.00001f && a - d > -0.00001f) {
      v1x = v2y = 1.0;
      v1y = v2x = 0.0;
      l1 = l2 = a;
    } else {
      v1x = -b/(a-d);
      v1y = 1.0;
      v2x = 1.0;
      v2y = 0.0;
      l1 = d;
      l2 = a;
    }
  } else {
    double ad = (a + d)/2.0;
    double s = sqrt(a*a - 2.0*a*d + d*d + 4.0*b*c)/2.0;
    v1x = (ad + s - d)/c;
    v1y = 1.0;
    v2x = (ad - s - d)/c;
    v2y = 1.0;
    l1 = ad + s;
    l2 = ad - s;
  }
  // Normalize
  double n = sqrtf(v1x*v1x + v1y*v1y);
  v1x /= n; v1y /= n;
  n = sqrtf(v2x*v2x + v2y*v2y);
  v2x /= n; v2y /= n;

  V->v[0] = v1x; V->v[1] = v2x;
  V->v[2] = v1y; V->v[3] = v2y;
  D->v[0] = l1;
  D->v[1] = l2;
}

void matrix_print(matrix_t *a, char *name) {
  int i, j;
  int m = a->m;
  int n = a->n;

  printf("%s =\t| ", name);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      //printf("%10.4f ", a->v[i*n+j]);
      printf("%13e ", a->v[i*n+j]);
    }
    printf("|\n");
    if (i != m-1) printf("\t| ");
  }
}

void array_print(float *a, char *name, int m, int n) {
  int i, j;

  printf("%s =\t| ", name);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("%10.4f ", a[i*n+j]);
    }
    printf("|\n");
    if (i != m-1) printf("\t| ");
  }
}

void transpose_matrix_print(matrix_t *a, char *name) {
  int i, j;
  int m = a->m;
  int n = a->n;

  printf("%s =\t| ", name);
  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      printf("%8.2f ", a->v[j*n+i]);
    }
    printf("|\n");
    if (i != n-1) printf("\t| ");
  }
}
