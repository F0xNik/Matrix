#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int ret = OK;
  if (rows < 1 || columns < 1) {
    ret = MATRIX_EROR;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = calloc(rows, sizeof(double *));
    if (result->matrix != NULL) {
      for (int i = 0; i < rows; i++) {
        result->matrix[i] = calloc(columns, sizeof(double));
      }
    } else
      ret = MATRIX_EROR;
  }
  return ret;
}

void fill_matrix(matrix_t *A) {
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      printf("El %d,%d: ", i, j);
      scanf("%lf", &A->matrix[i][j]);
    }
  }
}

void print_matrix(matrix_t A) {
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < A.columns; j++) {
      printf("%lf\t", A.matrix[i][j]);
    }
    printf("\n");
  }
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix) {
    for (int i = 0; i < A->rows; i++) {
      if (A->matrix[i]) free(A->matrix[i]);
    }
    free(A->matrix);
    A->matrix = NULL;
  }
  if (A->rows) A->rows = 0;
  if (A->columns) A->columns = 0;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int ret = 1;
  if (!is_exist(*A) || !is_exist(*B))
    ret = 0;
  else if (!is_size_eq(*A, *B))
    ret = 0;
  else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1E-7) ret = 0;
      }
    }
  }
  return ret;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int ret = OK;
  if (!is_exist(*A) || !is_exist(*B))
    ret = MATRIX_EROR;
  else if (!is_size_eq(*A, *B))
    ret = CALC_EROR;
  else {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else
      ret = MATRIX_EROR;
  }
  return ret;
}
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int ret = OK;
  if (!is_exist(*A) || !is_exist(*B))
    ret = MATRIX_EROR;
  else if (!is_size_eq(*A, *B))
    ret = CALC_EROR;
  else {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    } else
      ret = MATRIX_EROR;
  }
  return ret;
}

int is_exist(matrix_t A) {
  int ret = 1;
  if (A.rows < 1 || A.columns < 1 || A.matrix == NULL) ret = 0;
  return ret;
}

int is_size_eq(matrix_t A, matrix_t B) {
  int ret = 1;
  if (A.rows != B.rows || A.columns != B.columns) ret = 0;
  return ret;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int ret = OK;
  if (!is_exist(*A))
    ret = MATRIX_EROR;
  else {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    } else
      ret = MATRIX_EROR;
  }
  return ret;
}
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int ret = OK;
  if (!is_exist(*A) || !is_exist(*B))
    ret = MATRIX_EROR;
  else if (A->columns != B->rows)
    ret = CALC_EROR;
  else {
    if (s21_create_matrix(A->rows, B->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          for (int k = 0; k < B->rows; k++)
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    } else
      ret = MATRIX_EROR;
  }
  return ret;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int ret = OK;
  if (!is_exist(*A))
    ret = MATRIX_EROR;
  else {
    if (s21_create_matrix(A->columns, A->rows, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[j][i] = A->matrix[i][j];
        }
      }
    } else
      ret = MATRIX_EROR;
  }
  return ret;
}

int s21_determinant(matrix_t *A, double *result) {
  int ret = OK;
  if (!is_exist(*A))
    ret = MATRIX_EROR;
  else if (A->rows != A->columns)
    ret = CALC_EROR;
  else {
    if (A->rows == 1 && A->columns == 1)
      *result = A->matrix[0][0];
    else if (A->rows == 2 && A->columns == 2)
      *result = ((A->matrix[0][0] * A->matrix[1][1]) -
                 (A->matrix[1][0] * A->matrix[0][1]));
    else {
      double res = 0;
      for (int i = 0; i < A->rows; i++) {
        matrix_t minor;
        Minor(i, 0, A, &minor);
        s21_determinant(&minor, &res);
        *result += pow(-1, i) * A->matrix[0][i] * res;
        s21_remove_matrix(&minor);
        res = 0;
      }
    }
  }
  return ret;
}

void Minor(int i, int j, matrix_t *A, matrix_t *minor) {
  int k_m = 0;
  int l_m = 0;
  s21_create_matrix(A->rows - 1, A->columns - 1, minor);
  for (int k = 0; k < A->rows; k++) {
    if (k != j) {
      for (int l = 0; l < A->columns; l++) {
        if (l != i) {
          minor->matrix[k_m][l_m] = A->matrix[k][l];
          l_m++;
        }
      }
      l_m = 0;
      k_m++;
    }
  }
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int ret = OK;
  if (!is_exist(*A)) {
    ret = MATRIX_EROR;
  } else if (A->rows != A->columns) {
    ret = CALC_EROR;
  } else if (A->rows == 1 && A->columns == 1) {
    s21_create_matrix(A->rows, A->columns, result);
    result->matrix[0][0] = A->matrix[0][0];
  } else {
    if (A->columns != 1 && A->rows != 1) {
      double res = 0;
      matrix_t res_copy;
      s21_create_matrix(A->rows, A->columns, &res_copy);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          matrix_t minor;
          Minor(i, j, A, &minor);
          s21_determinant(&minor, &res);
          s21_remove_matrix(&minor);
          res_copy.matrix[i][j] = pow(-1, i + j) * res;
          res = 0;
        }
      }
      s21_transpose(&res_copy, result);
      s21_remove_matrix(&res_copy);
    } else {
      result->matrix[0][0] = A->matrix[0][0];
    }
  }
  return ret;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int ret = OK;
  if (!is_exist(*A))
    ret = MATRIX_EROR;
  else if (A->rows != A->columns)
    ret = CALC_EROR;
  else {
    double det = 0;
    s21_determinant(A, &det);
    matrix_t memory, kek;
    s21_calc_complements(A, &memory);
    s21_transpose(&memory, &kek);
    s21_mult_number(&kek, 1.0 / det, result);
    s21_remove_matrix(&memory);
    s21_remove_matrix(&kek);
  }
  return ret;
}