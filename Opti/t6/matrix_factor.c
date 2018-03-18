#include <math.h>
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