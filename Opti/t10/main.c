#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "optimization.h"

double* function(double *x, int n, double *f){
  if(n % 2) exit(1);
  for (int i = 0; 2 * i < n; i += 1) {
    f[2 * i] = 10 * (x[2 * i + 1] - SQUARE(x[2 * i]));
    f[2 * i + 1]   = 1 - x[2*i];
  }
  return f;
}
double **jacobian(double *x, int n, int m, double **j){
  //assume matrix comes all zeros...
  for(int i = 0; 2 * i < n; i++) {
    j[2 * i][2 * i]     = -20.0 * x[2 * i];
    j[2 * i + 1][2 * i] = -1.0;
    j[2 * i][2 * i + 1] = 10.0;
  }
  return j;
}
int main(int argc, char **argv){
  int n;
  FILE *file = fopen(argv[1], "r");
  fscanf(file, "%d", &n);
  double *x = newVector(n);
  for (int i = 0; i < n; ++i) {
    fscanf(file, "%lg", &x[i]);
  }
  fclose(file);
  int max_iter = atoi(argv[2]);
  double toler = atof(argv[3]);
  double v     = atof(argv[4]);
  // printf("%g\n", v);
  double *x_out = levenbergMarquardt(x, function, jacobian, max_iter, toler, n, n, v);
  // printVector(x_out, n);
  free(x_out);
  free(x);
}