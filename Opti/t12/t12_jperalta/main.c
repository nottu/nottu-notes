#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "optimization.h"

double f_eval(double *x, int n){
  return SQUARE(x[0] - 6) + SQUARE(x[1] - 7);
}
double h_1(double *x, int n){
  return -3 * x[0] - 2 * x[1] + 6;
}
double h_2(double *x, int n){
  return -x[0] + x[1] - 3;
}
double h_3(double *x, int n){
  return x[0] + x[1] - 7;
}
double h_4(double *x, int n){
  return (2.0/3.0) * x[0] - x[1] - (4.0/3.0);
}

int main(int argc, char **argv){
  double *x = newVector(2);
  switch(atoi(argv[1])){
    case 0:
      x[0] = -15; x[1] = 25;
      break;
    default:
      x[0] = -15; x[1] = -10;
  }
  FunEval *h = (FunEval*)malloc(sizeof(FunEval) * 4);
  h[0] = h_1; h[1] = h_2; h[2] = h_3; h[3] = h_4;
  double *x2 = quadraticPenalizaition(1, x, 0.005, 1000, 2, f_eval, NULL, h, 0, 4);
  printf("Aaaaaa\n");
  printf("%d\tx : (", i);printVector(x2, n);printf(")\t f:%g\t", -f(x2, n));
  free(x2);
  free(x);
}