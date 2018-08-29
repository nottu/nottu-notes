//
// Created by Javier Peralta on 9/16/17.
//

#ifndef TAREA3_MATRIX_H
#define TAREA3_MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MTXMAXSIZE 10000000 //~10MB
#define ZERO 0.000001

void vectorScalar(double *v, double d, int size);
void  restaVector(double *v1, double *v2, double* out, int size);
void   sumaVector(double *v1, double *v2, double* out, int size);
void   multVector(double *v1, double *v2, double* out, int size);
double sumElemetsVec(double *v1, int size);
//
void   productoPuntoA(double *v1, double *v2, double* out, int size);
double  productoPunto(double *v1, double *v2, int size);

void  sumaMatrix(double **mat1, double **mat2, int n, int m, double** out);
void restaMatrix(double **mat1, double **mat2, int n, int m, double** out);
//
void transposeMtx (double **mat, double** mt, int n, int m);
void multMatrizVect(double **mat,  double *vec,   int n, int m, double* res);
void     multMatriz(double **mat1, double **mat2, int n, int m, int p, int q, double** out);
void  multMatrizTrans(double **mat1, int n, double** res);
//
void printVect(double *a,  int n);
void  printMtx(double **a, int nr, int nc);
void  printMtxT(double **a, int nr, int nc);
//
double*  readVector(char* name, int* sz);
double**    readMtx(char* name, int* nr, int* nc);
//
double**  allocMtx(int nr, int nc);
void cleanMtx(double** mtx, int n, int m);
double**  allocMtxClean(int nr, int nc);
double** allocMtxI(int n);
void      freeMtx(double**a);
//
double   norma2Vect(double *v,   int size);
double norma2VectSq(double* v, int size);
void  normalizaVect(double *v, int size);
//
double   diffVectSq(double *v1,  double *v2,  int size);
double diffMatrizSq(double **m1, double **m2, int nr, int nc);
double     normaInf(double **m1, int n, int m);

//sol mtx
double*  diagSol(double  *a, double *b,  int n);
double* upperSol(double **a, double *b,  int nr, int nc);
double* lowerSol(double **a, double *b,  int nr, int nc);
//
double** cholesky(double **a, int n, int m);
void makedefpos(double **mat, int n);
int isdefpos(double **mat, int n);
int      luFactor(double **a, double **l, double **u, int nr, int nc);
double*  luSolver(double **l, double **u, double *b, int nr, int nc);
int     luFactor2(double **a, int nr, int nc);
double* luSolver2(double **a, double *b, int nr, int nc);

//sol tridiag
double* triDiagSol(double **a, double *b, int nr);
void printTridiag(double **d, int sz);
//eigen
double     potencia(double **mat, double *eigvec, int nr, int nc, int maxIter, double toler);
double  potenciaInv(double **mat, double *eigvec, double val, int nr, int nc, int maxIter, double toler, int *k, double *err);
double smallestEigv(double **mat, double *eigvec, int n, int m, int maxIter, double toler);
double  nearestEigv(double **mat, double *eigvec, double val,  int n, int m, int maxIter, double toler);
double*     allEigv(double **mat, int n, int m, int maxIter, double toler, int sections);
//inverse, mat is nxm, inv is mxn
void inverseMtx(double **mat, double **inv, int n, int m);

//jacobi
double valMayor(double **mat, int n, int m, int *x, int *y);
double* jacobiEig(double **mat,  double**eigVec, int n, int m, int maxIter, double toler);

//gradiente conjugado
//square matrix only
void gradienteConjugado(double **A, double *b, double *x, int n, double t);

//Polinomial interpolation

double** genVandermondeMtx(double *xvals, int n);
double evalPolynomial(const double *c, double x, int n);
double evalLagrange(double x, const double *xvals, int i, int n);
double evalPolynomialLagrange(double x, double *xvals, const double *yvals, int n);
double newtonPolynomial(double x, const double *xvals, int i);
double** newtonPolynomialMtx(double *xvals, const double *yvals, int n);
double evalNewtonPolynomial(double x, double *a, double *xvals, const double *yvals, int n);
#endif //TAREA3_MATRIX_H
