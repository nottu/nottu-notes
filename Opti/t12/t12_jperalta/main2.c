#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "optimization.h"

double f_eval(double *x, int n){
  return - pow(x[0], 2.0/3.0)*pow(x[1], 1.0/3.0);
}
double g_1(double *x, int n){
  return x[0] + x[1] - 4000;
}

int main(int argc, char **argv){
  int n = 2;
  double *x = newVector(n);
  x[0] = 40; x[1] = 100;
  FunEval *g = (FunEval*)malloc(sizeof(FunEval) * 1);
  g[0] = g_1;
  double *x2 = quadraticPenalizaition(1, x, 0.005, 1000, 2, f_eval, g, NULL, 1, 0);
  printf("x : (");printVector(x2, n);printf(")\t f:%g\t", -f_eval(x2, n));

  free(x2);
  free(x);
}