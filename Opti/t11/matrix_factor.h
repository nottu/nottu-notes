#ifndef MATRIX_FACTOR_H
#define MATRIX_FACTOR_H

int cholesky(double **a, int n, double **l);
int luFactor(double** a, double **l, double **u, int nr, int nc);

//tries cholesky, 1 if true 0 if false
int isdefpos(double **mat, int n);
void makedefpos(double **mat, int n);
double* luSolver(double **l, double **u, double *b, int nr, int nc);
double* luSolver2(double **a, double *b, int nr, int nc);
double matrix_condition(double **mat, int n);

//inverse, mat is nxm, inv is mxn
void inverseMtx(double **mat, double **inv, int n, int m);
#endif