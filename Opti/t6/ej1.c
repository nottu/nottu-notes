#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "optimization.h"

int doDoglegRoss(char* filename, int max_iter, double tg, double tx, double tf, double reg_sz){
  FuncInfo f;
  f.function = rosembrock;
  f.gradient = rosembrock_gradient;
  f.hessian  = rosembrock_hessian;

  int n;

  FILE *file = fopen(filename, "r");
  if(file == NULL) return 1; //could not open!
  fscanf(file, "%d", &n);
  double *x = newVector(n);
  for (int i = 0; i < n; ++i)
  {
    fscanf(file, "%lg", &x[i]);
  }
  fclose(file);
  doglegOptimize(f, x, n, max_iter, tg, tx, tf, reg_sz);
  free(x);

  return 0;
}

int main(int argc, char **argv){

  if(argc < 9){
    printf("t5 [rosenbrock] [Dogleg] [filename] [maxiter] [tol_gradient] [tol_x] [tol_function] [region_size]\n");
    return 1;
  }

  char *func_name = argv[1];
  char *step_type = argv[2]; //ignore for now
  char *file_name = argv[3];

  int max_iter = atoi(argv[4]);
  double tol_g = atof(argv[5]);
  double tol_x = atof(argv[6]);
  double tol_f = atof(argv[7]);

  double reg_sz = atof(argv[8]);

  return doDoglegRoss(file_name, max_iter, tol_g, tol_x, tol_f, reg_sz);
}