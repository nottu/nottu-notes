#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "optimization.h"


int main(int argc, char **argv){
  FuncInfo f;
  f.function = rosembrock;
  f.gradient = rosembrock_gradient;

  int n;
  FILE *file = fopen(argv[1], "r");
  fscanf(file, "%d", &n);
  double *x = newVector(n);
  for (int i = 0; i < n; ++i)
  {
    fscanf(file, "%lg", &x[i]);
  }
  fclose(file);

  // double **A = rosembrock_hessian(x, n, allocMatrix(n, n));
  double **A = HessianAproximator(x, f, n, atof(argv[2]));
  // printMatrix(A, n, n);
  double **H = allocMatrix(n, n);
  inverseMtx(A, H, n, n);
  // printMatrix(H, n, n);

  double *sol = bfgs(x, f, H, atoi(argv[3]), sqrt(DBL_EPSILON), n);
  printVector(sol, n);

  free(sol);
  freeMatrix(H);
  freeMatrix(A);
  free(x);
}