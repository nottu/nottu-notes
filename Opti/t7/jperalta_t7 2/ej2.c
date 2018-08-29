#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "optimization.h"

void test_rosembrock(char* filename, double tg, double tx, double tf, BetaFun betaFun){
  FuncInfo f;
  f.function = rosembrock;
  f.gradient = rosembrock_gradient;
  f.hessian  = rosembrock_hessian;
  int n;
  FILE *file = fopen(filename, "r");
  if(file == NULL){
    exit(1);
  }
  fscanf(file, "%d", &n);
  double *x = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; ++i)
  {
    fscanf(file, "%lg", &x[i]);
  }
  fclose(file);
  non_linear_conjugateGradient(f, x, n, tg, tx, tf, betaFun);
  free(x);
}

void test_wood(char* filename, double tg, double tx, double tf, BetaFun betaFun){
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
  for (int i = 0; i < n; ++i){
    fscanf(file, "%lg", &x[i]);
  }
  fclose(file);
  non_linear_conjugateGradient(f, x, n, tg, tx, tf, betaFun);
  free(x);
}
void test_convex1(char* filename, double tg, double tx, double tf, BetaFun betaFun){
  FuncInfo f;
  f.function = convex1;
  f.gradient = convex1_gradient;
  f.hessian  = convex1_hessian;
  int n;
  FILE *file = fopen(filename, "r");
  if(file == NULL){
    exit(1);
  }
  fscanf(file, "%d", &n);
  double *x = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; ++i){
    fscanf(file, "%lg", &x[i]);
  }
  fclose(file);
  non_linear_conjugateGradient(f, x, n, tg, tx, tf, betaFun);
  free(x);
}
void test_convex2(char* filename, double tg, double tx, double tf, BetaFun betaFun){
  FuncInfo f;
  f.function = convex2;
  f.gradient = convex2_gradient;
  f.hessian  = convex2_hessian;
  int n;
  FILE *file = fopen(filename, "r");
  if(file == NULL){
    exit(1);
  }
  fscanf(file, "%d", &n);
  double *x = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; ++i){
    fscanf(file, "%lg", &x[i]);
  }
  fclose(file);
  non_linear_conjugateGradient(f, x, n, tg, tx, tf, betaFun);
  free(x);
}


int main(int argc, char **argv){

  if(argc < 7){
    printf("t7 [rosembrock | wood | convex1 | convex2] [fletcher | polak | hestenes] [filename] [tol_gradient] [tol_x] [tol_function]\n");
    return 1;
  }
  char *func_name = argv[1];
  char *method = argv[2];
  char *file_name = argv[3];

  double tol_g = atof(argv[4]);
  double tol_x = atof(argv[5]);
  double tol_f = atof(argv[6]);

  BetaFun betaFun = fletcher_reeves;

  if(strcmp(method, "hestenes") == 0){
    betaFun = hestenes_stiefel;
  } else if(strcmp(method, "polak") == 0){
    betaFun = polak_ribiere;
  } else if(strcmp(method, "fletcher")){
    betaFun = fletcher_reeves;
  }

  if(strcmp(func_name, "rosembrock") == 0){
    test_rosembrock(file_name, tol_g, tol_x, tol_f, betaFun);
  }
  else if(strcmp(func_name, "wood") == 0){
    test_wood(file_name, tol_g, tol_x, tol_f, betaFun);
  }
  else if(strcmp(func_name, "convex1") == 0){
    test_convex1(file_name, tol_g, tol_x, tol_f, betaFun);
  }
  else if(strcmp(func_name, "convex2") == 0){
    test_convex2(file_name, tol_g, tol_x, tol_f, betaFun);
  }
  return 0;
}