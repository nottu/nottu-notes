#ifndef MATRIX_FACTOR_H
#define MATRIX_FACTOR_H

int cholesky(double **a, int n, double **l);

//tries cholesky, 1 if true 0 if false
int isdefpos(double **mat, int n);
#endif