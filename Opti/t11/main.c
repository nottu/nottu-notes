#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "optimization.h"

double eval_fun(double *p, double x){
  return p[0] * x + p[1] + p[2] * exp(p[3] * pow(x - p[4], 2));
}

int main(int argc, char **argv){
  int n, m;
  double **data = readMtx(argv[1], &n, &m);
  double **dataT = transposeMatrix(data, n, m);

  double *x = dataT[0];
  double *y = dataT[1];
  m = 5;

  //Nested Functions, only in GCC, x and y should stay in scope....
  double* function(double *p, int n, double *f){
    for (int i = 0; i < n; i += 1) {
      f[i] = p[0] * x[i] + p[1] + p[2] * exp(p[3] * pow(x[i] - p[4], 2)) - y[i];
    }
    return f;
  }
  double **jacobian(double *p, int n, int m, double **j){
    for(int i = 0; i < m; i++) {
      j[i][0] = x[i];
      j[i][1] = 1.0;
      double tmp = pow(x[i] - p[4], 2);
      j[i][2] = exp(p[3] * tmp);
      j[i][3] = p[2] * exp(p[3] * tmp) * tmp;
      j[i][4] = -2.0 * p[2] * p[3] * exp(p[3] * tmp) * (x[i] - p[4]);
    }
    return j;
  }

  double h = atof(argv[4]);

  double **jacobian_approx(double *p, int n, int m, double **J){
    double *p_t = newVector(n);
    copyVector(p, p_t, n);
    double **J1 = allocMatrix(m, n);
    double *r = newVector(m);
    for (int i = 0; i < n; ++i){
      p_t[i] += + h;
      function(p_t, m, r);
      for (int j = 0; j < m; ++j){
        J1[j][i] = r[j];
      }
      p_t[i] -= h;
    }
    function(p, m, r);
    for (int i = 0; i < m; ++i){
      for (int j = 0; j < n; ++j){
        J[i][j] = (J1[i][j] - r[i]) / h;
      }
    }
    free(r);
    free(p_t);
    freeMatrix(J1);
    return J;
  }

  double *p = newVector(m);
  p[0] = p[1] = 0.0;
  p[2] =  15.0;
  p[3] = -2.0;
  p[4] =  1.0;

  double toler = 0.001*n;
  int max_iter = atoi(argv[2]);
  double v     = atof(argv[3]);

  double *p_out = levenbergMarquardt(p, function, jacobian_approx, max_iter, toler, m, n, v);
  printf("\n");printVector(p_out, m);

  free(p_out);
  free(p);
  freeMatrix(dataT);
  freeMatrix(data);
  return 0;
}