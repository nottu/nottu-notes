#ifndef MATRIX_FACTOR_H
#define MATRIX_FACTOR_H

int cholesky(double **a, int n, double **l);
int luFactor(double** a, double **l, double **u, int nr, int nc);

//tries cholesky, 1 if true 0 if false
int isdefpos(double **mat, int n);

double* luSolver(double **a, double *b, int nr, int nc);
#endif