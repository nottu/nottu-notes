#include <stdio.h>
#include <stdlib.h>

#include "optimization.h"

int main(int argc, char **argv){
  FuncInfo f;
  f.function = smoothing;
  f.gradient = smoothing_gradient;
  f.hessian  = smoothing_hessian;
  int n;

  FILE *file = fopen(argv[1], "r");
  double tol_g = atof(argv[2]);
  double tol_x = atof(argv[3]);
  double tol_f = atof(argv[4]);
  double lambda = atof(argv[5]);

  fscanf(file, "%d", &n);
  //we read y, x into 'x' vector
  double *x = (double*)malloc(sizeof(double) * (2*n+ 1));
  for (int i = 0; i < n; ++i)
  {
    fscanf(file, "%lg", &x[n+i]);//y comes first
  }
  for (int i = 0; i < n; ++i)
  {
    fscanf(file, "%lg", &x[i]);//y comes first
  }
  x[2*n] = lambda;
  fclose(file);
  linealConjugateGradient(f, x, n, tol_g, tol_x, tol_f);
}