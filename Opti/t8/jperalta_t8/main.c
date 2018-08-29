#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "optimization.h"

double* function(double *x, int n, double *f){
  f[0] = 3 * x[0] - cos(x[1] * x[2]) - 0.5;
  f[1] = SQUARE(x[0]) - 81 * SQUARE(x[1] + 0.1) + sin(x[2]) + 1.06;
  f[2] = exp(-x[0] * x[1]) + 20 * x[2] + (10 * M_PI -3 ) / 3;
  return f;
}

double **jacobian(double *x, int n, double **j){
  double t1 = -x[0] * x[1];
  double t2 = x[1] * x[2];

  j[0][0] = 3.0;
  j[0][1] = x[2] * sin(t2);
  j[0][2] = x[1] * sin(t2);
  j[1][0] = 2.0 * x[0];
  j[1][1] = -162.0 * (x[1] + 0.1);
  j[1][2] = cos(x[2]);
  j[2][0] = -x[1] * exp(t1);
  j[2][1] = -x[0] * exp(t1);
  j[2][2] = 20.0;

  return j;
}

double* x_0_cases(int c){
  double *x = newVector(3);
  switch(c){
    case 0:
      x[0] = x[1] = x[2] = 0.0;
      break;
    case 1:
      x[0] = x[1] = 1.1;
      x[2] = -1.1;
      break;
    case 2:
      x[0] = x[1] = -10.0;
      x[2] = 10.0;
      break;
    case 3:
      x[0] = x[2] = 3.0;
      x[1] = -3.0;
      break;
    default:
      x[0] = x[1] = x[2] = 0.0;
      break;
  }
  return x;
}


int main(int argc, char **argv){
  int met = atoi(argv[1]);
  double *x = x_0_cases( atoi(argv[3]) );
  int n = 3;

  int max_iter = atoi(argv[2]);
  double *x_out;
  if (met == 0)x_out = newton_no_lineal(x, function, jacobian, max_iter, sqrt(DBL_EPSILON), n);
  else x_out = broyden(x, function, jacobian, max_iter, sqrt(DBL_EPSILON), n);
  free(x_out);
  free(x);
}