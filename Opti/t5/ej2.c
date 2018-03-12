#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "optimization.h"

double number_pii(double *b, double *x, int n){
  double pi_i = 1.0 / (1.0 + exp(-dotproduct(x, b, n) - b[n]));
  return pi_i;
}

double number(double *b, int *y, double **x, int n, int m){
  double ev = 0;
  for (int i = 0; i < n; ++i){
    double pi_i = number_pii(b, x[i], m);
    if (y[i] == 0) ev += log(1 - pi_i);
    else ev += log(pi_i);
  }
  return -ev;
}

double *number_gradient(double *b, int *y, double **x, int n, int m, double *g){
  double *pi_is = newVector(n);
  for (int i = 0; i < n; ++i) {
    pi_is[i] = number_pii(b, x[i], m);
  }
  for (int j = 0; j < m; ++j) {
    g[j] = 0;
    for (int i = 0; i < n; ++i) {
      g[j] += x[i][j] * ( y[i] - pi_is[i]);
    }
  }
  g[m] = 0;
  for (int i = 0; i < n; ++i) {
    g[m] += y[i] - pi_is[i];
  }
  scaleVector(g, m+1, -1);
  free(pi_is);
  return g;
}

double number_step_backtracking(double* d, double *b, int *y, double **x, double c, double p, double alpha, int n, int m){
  double *g = newVector(m + 1);
  number_gradient(b, y, x, n, m, g);
  double fd = dotproduct(g, d, m + 1);
  double fx = number(b, y, x, n, m);
  int maxiter = 1000;
  while(maxiter--){
    memcpy(g, d, sizeof(double) * (m + 1)); //copy direction to gradient vector
    scaleVector(g, m+1, alpha);
    sumVectors(b, g, g, m+1);
    double fxk = number(g, y, x, n, m);
    if(fxk <= fx + alpha*c*fd) break;
    alpha *= p;
  }
  free(g);
  return maxiter ? alpha : 0;
}


void train_number(char* x_filename, char* y_filename, int num1, int num2){
  printf("TRAINING!\n");
  int num_train = 1000;
  int n_x = 28*28;
  FILE *yf = fopen(y_filename, "r"), *xf = fopen(x_filename, "r");
  int y_n = -1;
  double x_n = -1;

  int *y = (int*)malloc(sizeof(int) * num_train);
  double **x = allocMatrix(num_train, n_x);
  int cnt = 0;
  while(!feof(yf) && !feof(xf) && cnt < num_train){
    fscanf(yf, "%d", &y_n);
    if(y_n == num1 || y_n == num2){
      if(y_n == num1) y[cnt] = 0;
      if(y_n == num2) y[cnt] = 1;
      for (int i = 0; i < n_x -1; ++i){
        fscanf(xf, "%lf,", &x_n);
        x[cnt][i] = x_n;
      }
      fscanf(xf, "%lf", &x_n);
      x[cnt][n_x - 1] = x_n;
      cnt++;
    } else {
        for (int i = 0; i < n_x -1; ++i){
        fscanf(xf, "%lf,", &x_n);
      }
      fscanf(xf, "%lf", &x_n);
    }
  }
  fclose(yf);
  fclose(xf);

  double *beta = newVector(n_x + 1);
  double *beta1 = newVector(n_x + 1);

  for (int i = 0; i < n_x+1; ++i) beta[i] = 1.0/(n_x + 1);

  double ev = number(beta, y, x, num_train, n_x);
  printf("val inicial: %g\n", -ev);

  double *g = newVector(n_x+1);
  int i = 0;
  for (i = 0; i < 1000000; ++i)
  {
    double pev = ev;
    g = number_gradient(beta, y, x, num_train, n_x, g);
    // printVector(g, n_x+1);
    scaleVector(g, n_x+1, -1);
    double step = 0.01;
    step = number_step_backtracking(g, beta, y, x, 1E-4, 0.99, 1, num_train, n_x);
    if(step == 0) break;
    scaleVector(g, n_x+1, step);
    sumVectors(beta, g, beta1, n_x+1);
    ev = number(beta1, y, x, num_train, n_x);
    printf("EV: %g\t step: %g\n", ev, step);
    for (int j = 0; j < n_x+1; ++j) beta[j] = beta1[j];
    if(fabs(pev - ev) < 1E-8) break;
  }
  // printVector(beta1, n_x+1);
  FILE *bout = fopen("b.out", "w");
  for (int i = 0; i < n_x+1; ++i)
  {
    fprintf(bout, "%lf\n", beta[i]);
  }
  printf("#Iter: %i\nval final: %g\n",i, -ev);
  free(beta);
  freeMatrix(x);
  free(y);
}

void test_number(char* x_filename, char* y_filename, int num1, int num2){
  printf("TESTING!\n");
  int num_train = 5000;
  int n_x = 28*28;
  FILE *yf = fopen(y_filename, "r"), *xf = fopen(x_filename, "r");
  int y_n = -1;
  double x_n = -1;

  int *y = (int*)malloc(sizeof(int) * num_train);
  double **x = allocMatrix(num_train, n_x);
  int cnt = 0;
  while(!feof(yf) && !feof(xf) && cnt < num_train){
    fscanf(yf, "%d", &y_n);
    if(y_n == num1 || y_n == num2){
      if(y_n == num1) y[cnt] = 0;
      if(y_n == num2) y[cnt] = 1;
      for (int i = 0; i < n_x -1; ++i){
        fscanf(xf, "%lf,", &x_n);
        x[cnt][i] = x_n;
      }
      fscanf(xf, "%lf", &x_n);
      x[cnt][n_x - 1] = x_n;
      cnt++;
    } else {
        for (int i = 0; i < n_x -1; ++i){
        fscanf(xf, "%lf,", &x_n);
      }
      fscanf(xf, "%lf", &x_n);
    }
  }
  fclose(yf);
  fclose(xf);

  FILE *b_file = fopen("b.out", "r");
  double *b = newVector(n_x+1);
  for (int i = 0; i < n_x+1; ++i){
    fscanf(b_file, "%lf,", &b[i]);
  }

  double err = 0;
  for (int i = 0; i < num_train; ++i){
    double pi_i = number_pii(b, x[i], n_x+1);
    if((pi_i < 0.5 && y[i] == 1) || (pi_i > 0.5 && y[i] == 0)) err += 1.0;
  }
  err /= num_train;
  printf("ERROR %lg\n", err);
}

int main(int argc, char **argv){

  if (argc < 6)
  {
    printf("ej2 [trainY] [trainX] [testY] [testX] [num1] [num2]\n");
    return 1;
  }

  int num1 = atoi(argv[5]);
  int num2 = atoi(argv[6]);

  train_number(argv[2], argv[1], num1, num2);
  test_number(argv[4], argv[3], num1, num2);
}