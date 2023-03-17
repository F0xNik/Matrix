#ifndef _S21_MATRIX_H_
#define _S21_MATRIX_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SUCCESS 1
#define FAILURE 0

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

enum func_state { OK = 0, MATRIX_EROR = 1, CALC_EROR = 2 };

int s21_create_matrix(int rows, int columns, matrix_t *result);
void print_matrix(matrix_t result);
void fill_matrix(matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int is_exist(matrix_t A);
int is_size_eq(matrix_t A, matrix_t B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
void Minor(int i, int j, matrix_t *A, matrix_t *minor);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

#endif  //  _S21_MATRIX_H