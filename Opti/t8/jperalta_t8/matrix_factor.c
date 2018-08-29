#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_factor.h"
#include "matrix.h"

#define ZEROISH 1E-6

int cholesky(double **a, int n, double **l){
  for (int i = 0; i < n; ++i) {
    l[i][i] = a[i][i];
    for (int k = 0; k < i; ++k) {
      l[i][i] -= l[i][k] * l[i][k];
    }
    if(l[i][i] <= ZEROISH){
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
      if(fabs(l[i][i]) < ZEROISH)
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
double* luSolver(double **l, double **u, double *b, int nr, int nc){
  double* sol = lowerSol(l, b, nr, nc);
  double* sol2 = upperSol(u, sol, nr, nc);
  free(sol);
  return sol2;
}
double* luSolver2(double **a, double *b, int nr, int nc){ //should move to another file!
  double **l = allocMatrix(nr, nc);
  double **u = allocMatrix(nr, nc);
  if(!luFactor(a, l, u, nr, nc)){
    printf("FAILED LU FACTOIRZATION!!\n");
    return NULL;
  }
  double* sol = lowerSol(l, b, nr, nc);
  double* sol2 = upperSol(u, sol, nr, nc);
  free(sol);
  freeMatrix(l);
  freeMatrix(u);
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
      if(fabs(a[i][i]) < ZEROISH)
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
void inverseMtx(double **mat, double **inv, int n, int m){
  double **l = allocMatrix(n, m);
  double **u = allocMatrix(n, m);
  if (luFactor(mat, l, u, n, m)){
    double *b = newVector(m);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < m; ++j) {
        b[j] = j == i;
      }
      double *sol = luSolver(l, u, b, n, m);
      for (int j = 0; j < n; ++j) {
        inv[j][i] = sol[j];
      }
      free(sol);
    }
    free(b);
  }
  freeMatrix(l);
  freeMatrix(u);
}

double max_value(double **mat, int n, int m, int *x, int *y){
  double mayor = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      if (i == j) continue;
      if(mayor < fabs(mat[i][j])){
        mayor = fabs(mat[i][j]);
        *x = i; *y = j;
      }
    }
  }
  return mayor;
}
//GT * A * G
void givensRotate(double **mat, int n, int m, int mi, int mj, double c, double s){
  for (int i = 0; i < m; ++i) {
    double matimi = mat[i][mi];
    mat[i][mi] = matimi * c - s*mat[i][mj];
    mat[i][mj] = matimi * s + c*mat[i][mj];
  }
  for (int i = 0; i < n; ++i) {
    double matmii =  mat[mi][i];
    mat[mi][i] = mat[mi][i] * c - s*mat[mj][i];
    mat[mj][i] = matmii * s + c*mat[mj][i];
  }
}
void givensM(double **mat, int n, int m, int mi, int mj, double c, double s){
  for (int i = 0; i < m; ++i) {
    double matimi = mat[i][mi];
    mat[i][mi] = matimi * c - s*mat[i][mj];
    mat[i][mj] = matimi * s + c*mat[i][mj];
  }
}

double* jacobiEig(double **mat, double**eigVec, int n, int m, int maxIter, double toler){
  int x, y;
  double max = max_value(mat, n, m, &x, &y);
  if(max < toler) return NULL; //eigvs in diag
  double **eigvalsM = allocMatrix(n, m);
  for (int i = 0; i < n; ++i) memcpy(eigvalsM[i], mat[i], sizeof(double) * m);
  int iter = 0;
  while (max > toler && ++iter < maxIter){
    double d = (eigvalsM[y][y] - eigvalsM[x][x])/(2 * eigvalsM[x][y]);
    double t = 1 / (fabs(d) + sqrt(1 + d*d));
    t = d > 0 ? t : -t;
    double c = 1/(sqrt(1 + t * t));
    double s = c * t;
    givensRotate(eigvalsM, n, m, x, y, c, s);
    givensM(eigVec, n, n, x, y, c, s);
    max = max_value(eigvalsM, n, m, &x, &y);
  }
  double **AV = allocMatrix(n, m);
  multiplyMatrix(mat, eigVec, n, n, n, AV);
  double **VD = allocMatrix(n, m);
  multiplyMatrix(eigVec, eigvalsM, n, n, n, VD);

  freeMatrix(VD);
  freeMatrix(AV);

  double *eigvals = newVector(n);
  for (int i = 0; i < n; ++i) {
    eigvals[i] = eigvalsM[i][i];
  }
  freeMatrix(eigvalsM);
  return eigvals;
}


double matrix_condition(double **mat, int n){
  double **eigMat = allocMatrixIdentity(n, n);
  double *eigvals = jacobiEig(mat, eigMat, n, n, 10000, 1E-8);
  // printVector(eigvals, n);
  freeMatrix(eigMat);
  double small = eigvals[0];
  double big = eigvals[0];
  for (int i = 1; i < n; ++i){
    if(fabs(eigvals[i]) < fabs(small)) small = fabs(eigvals[i]);
    else if(fabs(eigvals[i]) > fabs(big)) big = fabs(eigvals[i]);
  }
  free(eigvals);
  return big/small;
}