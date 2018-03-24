#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix_factor.h"
#include "matrix.h"

int cholesky(double **a, int n, double **l){
  for (int i = 0; i < n; ++i) {
    l[i][i] = a[i][i];
    for (int k = 0; k < i; ++k) {
      l[i][i] -= l[i][k] * l[i][k];
    }
    if(l[i][i] <= 1E-6){
      // printf("No se puede hacer factorización Cholesky\n");
      return 0; //could not do cholesky
    }
    l[i][i] = sqrt(l[i][i]);
    for (int j = i+1; j < n; ++j) {
      l[i][j] = 0;
      double lij = a[j][i];
      for (int k = 0; k < i; ++k) {
        lij -= l[j][k]*l[i][k];
      }
      lij /= l[i][i];
      l[j][i] = (lij);
    }
  }
  return 1;
}
//tries cholesky, 1 if true 0 if false
int isdefpos(double **mat, int n){
  double **l = allocMatrix(n, n);
  int defpos = cholesky(mat, n, l);
  freeMatrix(l);
  return defpos;
}
double* upperSol(double**a , double*b, int nr, int nc){
  double *vect = malloc(sizeof(double) * nc);
  for (int i = nr -1; i >= 0; i--) {
    double tmp = b[i];
    for (int j = i+1; j < nc; ++j) {
      tmp -= vect[j] * a[i][j];
    }
    vect[i] = tmp / a[i][i];
  }
  return vect;
}
double* lowerSol(double**a , double*b, int nr, int nc){
  double *vect = malloc(sizeof(double) * nc);

  for (int i = 0; i < nr; ++i) {
    double tmp = b[i];
    for (int j = 0; j < i && j < nc; ++j) {
      tmp -= vect[j] * a[i][j];
    }
    tmp /= a[i][i];
    vect[i] = tmp;
  }

  return vect;
}
int luFactor(double** a, double **l, double **u, int nr, int nc){
  for (int i = 0; i < nr; ++i) {
    u[i][i] = 1;
    for (int j = 0; j <= i && j <nc; ++j) {
      double lij = a[i][j];
      for (int k = 0; k < j; ++k) {
        lij -= l[i][k]*u[k][j];
      }
      l[i][j] = lij;
    }
    for (int j = i+1; j < nc; ++j) {
      double lij = a[i][j];
      if(fabs(l[i][i]) < 1E-4)
        return 0;
      for (int k = 0; k < i; ++k) {
        lij -= l[i][k]*u[k][j];
      }
      lij /= l[i][i];
      u[i][j] = lij;
    }
  }
  return 1;
}
double* luSolver(double **a, double *b, int nr, int nc){ //should move to another file!
  double **l = allocMatrix(nr, nc);
  double **u = allocMatrix(nr, nc);
  if(!luFactor(a, l, u, nr, nc)){
    printf("FAILED LU FACTOIRZATION!!\n");
    exit(1);
  }
  double* sol = lowerSol(l, b, nr, nc);
  double* sol2 = upperSol(u, sol, nr, nc);
  free(sol);
  return sol2;
}
//same as lu factor, but in 1 matrix
int luFactor2(double **a, int nr, int nc){
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j <= i && j <nc; ++j) {
      double lij = a[i][j];
      for (int k = 0; k < j; ++k) {
        lij -= a[i][k]*a[k][j];
      }
      a[i][j] = lij;
    }
    for (int j = i+1; j < nc; ++j) {
      double lij = a[i][j];
      if(fabs(a[i][i]) < 1E-4)
        return 0;
      for (int k = 0; k < i; ++k) {
        lij -= a[i][k]*a[k][j];
      }
      lij /= a[i][i];
      a[i][j] = lij;
    }
  }
  return 1;
}