#ifndef MATRIX_H
#define MATRIX_H
#define MTXMAXSIZE 10000000 //~10MB
#include "matrix_factor.h"
#include "vector.h"
//memory management
double** allocMatrix(int nr, int nc);
double** allocMatrixClean(int nr, int nc);
double** allocMatrixIdentity(int nr, int nc);
void freeMatrix(double **mtx);

//simple operations
void addMatrix(double **m1, double** m2, int n, int m, double **out);
void substractMatrix(double **m1, double** m2, int n, int m, double **out);
//multiply nxm with mxp
void multiplyMatrix(double **m1, double **m2, int n, int m, int p, double **out);
void multiplyMatrixTransposed(double **m1, double **m2, int n, int m, int p, double **out);
//vector matrix
void multiplyMatrixVector(double **mat, double *v, int n, int m, double *out);

void printMatrix(double**a, int nr, int nc);

#endif