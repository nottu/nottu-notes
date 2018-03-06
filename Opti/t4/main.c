#include <stdio.h>
#include <string.h>
#include <math.h>

#include "matrix.h"
#include "optimization.h"

void test_rosembrock(char* filename, int max_iter, double tg, double tx, double tf, Step step, double alpa_i){
  //read from file...
  FuncInfo f;
  f.function = rosembrock;
  f.gradient = rosembrock_gradient;
  f.hessian  = rosembrock_hessian;

  int n;

  FILE *file = fopen(filename, "r");
  fscanf(file, "%d", &n);
  double *x = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; ++i)
  {
    fscanf(file, "%lg", &x[i]);
  }
  fclose(file);
  optimize_function(f, step, x, n, max_iter, tg, tx, tf, alpa_i);
  free(x);
}


void test_wood(char* filename, int max_iter, double tg, double tx, double tf, Step step, double alpa_i){
  //read from file...
  FuncInfo f;
  f.function = wood;
  f.gradient = wood_gradient;
  f.hessian  = wood_hessian;

  int n;

  FILE *file = fopen(filename, "r");
  fscanf(file, "%d", &n);
  if(n != 4){
    printf("Invalid number of params for wood function\n");
    return;
  }
  double *x = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; ++i)
  {
    fscanf(file, "%lg", &x[i]);
  }
  fclose(file);
  optimize_function(f, step, x, n, max_iter, tg, tx, tf, alpa_i);
  free(x);
}

void test_smooth(char* filename, int max_iter, double tg, double tx, double tf, Step step, double alpa_i, double lambda){
  //read from file...
  FuncInfo f;
  f.function = smoothing;
  f.gradient = smoothing_gradient;
  f.hessian  = smoothing_hessian;

  int n;

  FILE *file = fopen(filename, "r");
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
  // printVect(x, n);
  optimize_function(f, step, x, n, max_iter, tg, tx, tf, alpa_i);
  free(x);
}


int main(int argc, char **argv){

  if(argc < 10){
    printf("t4 [rosenbrock | wood | smoothing] [StepFijo | StepHess | StepAprox] [filename] [maxiter] [tol_gradient] [tol_x] [tol_function] [initial_alpha] [lambda]\n");
    return 1;
  }
  char *func_name = argv[1];
  char *step_type = argv[2];
  char *file_name = argv[3];

  int max_iter = atoi(argv[4]);
  double tol_g = atof(argv[5]);
  double tol_x = atof(argv[6]);
  double tol_f = atof(argv[7]);

  double alp_i = atof(argv[8]);
  double lambda = atof(argv[9]);

  Step step = StepFijo;

  if(strcmp(step_type, "StepAprox") == 0){
    step = StepAprox;
  } else if(strcmp(step_type, "StepHess") == 0){
    step = StepHess;
  } else if(strcmp(step_type, "StepFijo")){
    printf("Tipo de paso %s no existe, se usara StepFijo\n", step_type);
  }

  if(strcmp(func_name, "rosenbrock") == 0){
    test_rosembrock(file_name, max_iter, tol_g, tol_x, tol_f, step, alp_i);
  }
  else if(strcmp(func_name, "wood") == 0){
    test_wood(file_name, max_iter, tol_g, tol_x, tol_f, step, alp_i);
  }
  else if(strcmp(func_name, "smoothing") == 0){
    test_smooth(file_name, max_iter, tol_g, tol_x, tol_f, step, alp_i, lambda);
  }


  // ros 2
  // n = 2;
  // f.function = rosembrock;
  // f.gradient = rosembrock_gradient;
  // f.hessian  = rosembrock_hessian;
  // double x[] = {-1.2, 1};

  //ros 100
  // n = 100;
  // f.function = rosembrock;
  // f.gradient = rosembrock_gradient;
  // f.hessian  = rosembrock_hessian;
  // double x[n];
  // for (int i = 0; i < n; ++i) x[i] = 1;
  // x[0] = -1.2;
  // x[98] = -1.2;

  //wood
  // n = 4;
  // f.function = wood;
  // f.gradient = wood_gradient;
  // f.hessian = wood_hessian;
  // double x[] = {-3, -1, -3, -1};

  // optimize_function(f, StepHess, x, n, 100000, 1e-6, 1e-10, 1e-12);
  // test_rosembrock();
  return 0;
}