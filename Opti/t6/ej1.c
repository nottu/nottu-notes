#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "optimization.h"

//almost the same as get_step_hess
double cauchy_point(FuncInfo info, double *x, int n, double *g, double **h){
  double gtg = dotproduct(g, g, n);

  info.hessian(x, n, h);
  double *hg = newVector(n);

  multiplyMatrixVector(h, g, n, n, hg);
  double alp = gtg / dotproduct(g, hg, n);

  free(hg);
  return -1 * alp;
}

void doglegDirection(FuncInfo info, double *x, int n, double *g, double **h, double* cauchy_dir, double reg_sz, double *dir){
  double *pb = luSolver(h, g, n, n);
  if(vectorNorm(pb, n, 2) < reg_sz) {
    // printf("SMALL DOG\n");
    copyVector(pb, dir, n);
    scaleVector(dir, n, -1);
    free(pb);
    return;
  }
  double *pdif = newVector(n);

  substractVectors(pb, cauchy_dir, pdif, n);

  double a = SQUARE(vectorNorm(pdif, n , 2)),
         b = 2 * dotproduct(pb, pdif, n), 
         c = SQUARE(vectorNorm(cauchy_dir, n, 2)) + SQUARE(reg_sz);
  double lamda = cuadraticSolver(a, b, c);

  scaleVector(pdif, n, lamda);
  sumVectors(cauchy_dir, pdif, dir, n);
  // copyVector(cauchy_dir, dir, n);
  free(pdif);
  free(pb);
}

double taylorEval(double fx, double *x, int n, double *g, double **h, double *p){
  double pg = dotproduct(g, p, n);
  double *hp = newVector(n);
  multiplyMatrixVector(h, p, n, n, hp);
  double php = dotproduct(hp, p, n);
  free(hp);

  return fx + pg + php/2.0;
}

double doglegOptimize(FuncInfo info, double *x, int n, int max_iter, double tg, double tx, double tf, double reg_szM){
  double *g  = info.gradient(x, n, newVector(n));
  double **H = info.hessian(x, n, allocMatrix(n, n));

  double val = info.function(x, n);
  double norm2 = vectorNorm(g, n, 2);
  double reg_sz = reg_szM;
  int iter = 0;

  double *pu = newVector(n);
  double *dir = newVector(n);
  double *x1 = newVector(n);
  while(norm2 >= tg && iter < max_iter){
    double fx = info.function(x, n);
    iter++;
    double a = cauchy_point(info, x, n, g, H);
    copyVector(g, pu, n);
    scaleVector(pu, n, a);
    double pu_norm = vectorNorm(pu, n ,2);
    copyVector(pu, dir, n);
    if(pu_norm < reg_sz){
      if(isdefpos(H, n)) {
        //do dogleg direction
        // printf("DOGLEG!\n");
        doglegDirection(info, x, n, g, H, pu, reg_sz, dir);
      }
    } else {
      //make smaller
      scaleVector(dir, n, reg_sz/pu_norm);
    }

    sumVectors(x, dir, x1, n); //x = x + dir
    //check if solution is good!
    double taylorev = fx - taylorEval(fx, x, n, g, H, dir);
    double phi = (fx - info.function(x1, n)) / (taylorev);
    printf("Iter %i:%i \tf(x): %g\t||g||: %g\t reg_sz %lg\t phi: %g\t taylor %g\n", iter, max_iter, fx, norm2, reg_sz, phi, taylorev);
    if(phi < 0.25) {
      reg_sz /= 4; //bad model
      // printf("\tBAD MODEL\tphi %g\n", phi);
    }
    if(phi < 0) continue;
    double dirnorm = vectorNorm(dir, n, 2);
    if(phi >= 0.75 && dirnorm + 0.001 > reg_sz ){
      printf("HERE!\t%g\t%g\n", reg_sz, 2*reg_sz);
      reg_sz = 2*reg_sz < reg_szM ? 2*reg_sz : reg_szM;
    }
    copyVector(x1, x, n);

    H = info.hessian(x, n, H);
    g = info.gradient(x, n, g);
    norm2 = vectorNorm(g, n, 2);
  }

  free(g);
  free(x1);
  freeMatrix(H);
  return val;
}

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