#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "optimization.h"

void test_rosembrock(char* filename, int max_iter, double tg, double tx, double tf, Step step, double alpa_i){
  FuncInfo f;
  f.function = rosembrock;
  f.gradient = rosembrock_gradient;
  f.hessian  = rosembrock_hessian;

  int n;

  FILE *file = fopen(filename, "r");
  fscanf(file, "%d", &n);
  double *x = newVector(n);
  for (int i = 0; i < n; ++i)
  {
    fscanf(file, "%lg", &x[i]);
  }
  fclose(file);
  optimize_function(f, step, x, n, max_iter, tg, tx, tf, alpa_i);
  free(x);
}


void test_wood(char* filename, int max_iter, double tg, double tx, double tf, Step step, double alpa_i){
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
  double *x = newVector(n);
  for (int i = 0; i < n; ++i)
  {
    fscanf(file, "%lg", &x[i]);
  }
  fclose(file);
  optimize_function(f, step, x, n, max_iter, tg, tx, tf, alpa_i);
  free(x);
}

int main(int argc, char **argv){

  if(argc < 9){
    printf("t5 [rosenbrock | wood | smoothing] [StepBacktrack | StepInterpol] [filename] [maxiter] [tol_gradient] [tol_x] [tol_function] [initial_alpha]\n");
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
  Step step = StepBacktrack;
  if(strcmp(step_type, "StepInterpol") == 0){
    step = StepInterpol;
  } else if(strcmp(step_type, "StepBacktrack")){
    printf("Tipo de paso %s no existe, se usara StepBacktrack\n", step_type);
  }

  if(strcmp(func_name, "rosenbrock") == 0){
    test_rosembrock(file_name, max_iter, tol_g, tol_x, tol_f, step, alp_i);
  }
  else if(strcmp(func_name, "wood") == 0){
    test_wood(file_name, max_iter, tol_g, tol_x, tol_f, step, alp_i);
  }
  return 0;
}