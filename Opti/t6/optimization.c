#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "optimization.h"

//step getting functions
double get_step_hess(FuncInfo info, double *x, int n, double *g){
  double gtg = dotproduct(g, g, n);
  double **h = allocMatrixClean(n, n);

  info.hessian(x, n, h);
  double *hg = newVector(n);
  
  multiplyMatrixVector(h, g, n, n, hg);
  double alp = gtg / dotproduct(g, hg, n);
  
  // if(!isdefpos(h, n)) alp *= -1;
  // makedefpos(h, n);
  free(hg);
  freeMatrix(h);
  return alp;
}
double get_step_approx(FuncInfo info, double *x, int n, double alp, double* g){
  double gtg = dotproduct(g, g, n);
  double numerator = gtg * alp * alp;

  double *dir = newVector(n);
  info.gradient(x, n, dir);
  double *x1  = newVector(n);
  scaleVector(dir, alp, n);
  substractVectors(x, dir, x1, n);

  double fx  = info.function(x1, n);
  double pfx = info.function(x,  n);

  double denominator = 2.0 * (fx - pfx + alp * gtg);
  free(dir);
  free(x1);
  return numerator/denominator;
}
double get_step_backtracking(FuncInfo info, double* d, double *x, double c, double p, double alpha, int n){
  double *g = newVector(n);
  info.gradient(x, n, g);
  double fd = dotproduct(g, d, n);
  double fx = info.function(x, n);
  while(1){
    memcpy(g, d, sizeof(double) * n); //copy direction to gradient vector
    scaleVector(g, n, alpha);
    sumVectors(x, g, g, n);
    double fxk = info.function(g, n);
    if(fxk <= fx + alpha*c*fd) break;
    alpha *= p;
  }
  free(g);
  return alpha;
}
//called by get_step_interpolated, not to be called by itself!!!
double __get_step_interpolated_cube(double a_0, double a_1, double phi_0, double phi_0p, double phi_a0, double phi_a1){
  double d = phi_0;
  double c = phi_0p;

  double div = SQUARE(a_0 * a_1) * (a_1 - a_0);
  double t1 = phi_a1 - phi_0p*a_1 - phi_0;
  double t2 = phi_a0 - phi_0p*a_0 - phi_0;

  double a = (t1 * SQUARE(a_0) - t2 * SQUARE(a_1))/div;
  double b = (-t1 * a_0*a_0*a_0 + t2 * a_1*a_1*a_1)/div;

  // printf("%g %g %g %g\n", a, b, c, d);

  b = -b + sqrt(SQUARE(b) - 3*a*c);

  return b/(3 * a);
}
double get_step_interpolated(FuncInfo info, double* d, double *x, double c, double alpha, int n){
  double *g = newVector(n);
  info.gradient(x, n, g);
  // printVector(g, n);
  double phi_0p = dotproduct(g, d, n);

  memcpy(g, d, sizeof(double) * n); //copy direction to gradient vector
  scaleVector(g, n, alpha);
  sumVectors(x, g, g, n);
  // printf("Alpha %g\n", alpha);

  double phi_0 = info.function(x, n);
  double phi_a = info.function(g, n);
  double alpha_1 = alpha;
  if(phi_a > phi_0 + c*alpha*phi_0p){
    //do interpolation stuff
    // printf("alpha :%g\tphi_0p: %g\tphi_a: %g\tphi_0: %g\n", alpha, phi_0p, phi_a, phi_0);
    alpha_1 = -SQUARE(alpha)*phi_0p/(2 * (phi_a - phi_0p*alpha - phi_0));
  }
  memcpy(g, d, sizeof(double) * n); //copy direction to gradient vector
  scaleVector(g, n, alpha);
  sumVectors(x, g, g, n);
  double phi_a1 = info.function(g, n);
  if(phi_a1 > phi_0 + c*alpha*phi_0p){
    // printf("Aun no cumple con Armijo!\t%g\n", alpha_1);
    alpha = __get_step_interpolated_cube(alpha, alpha_1, phi_0, phi_0p, phi_a, phi_a1);
    // printf("Alpha %g\n", alpha);
  } else {
    alpha = alpha_1;
  }

  free(g);
  return alpha;
}
double get_step(double alp, double alp_0, Step stp, FuncInfo info, double *gradient, double *x, int n){
  double step = alp;
  switch(stp){
    case StepFijo:
      break;
    case StepHess:
      step = get_step_hess(info, x, n, gradient);
      break;
    case StepAprox:
      step = get_step_approx(info, x, n, step, gradient);
      break;
    case StepBacktrack:
      scaleVector(gradient, n, -1);
      step = get_step_backtracking(info, gradient, x, 1E-4, 0.9, alp_0, n);
      scaleVector(gradient, n, -1);
      break;
    case StepInterpol:
      scaleVector(gradient, n, -1);
      step = get_step_interpolated(info, gradient, x, 1E-4, alp_0, n);
      scaleVector(gradient, n, -1);
    default:
      break;
  }
  return step;
}

//steepest decent
double optimize_function(FuncInfo info, Step stp, double *x, int n, int iter, double tg, double tx, double tf, double step){
  double *x1  = newVector(n);
  double *dir = newCleanVector(n);
  double fdif = 0;
  double fx = 0;
  FILE *file = fopen("data.txt", "w");
  double o_step = step;
  for (int i = 0; i < iter; ++i)
  {
    fx = info.function(x, n);
    info.gradient(x, n, dir); //max increse, use negative step...
    double norm_grad = vectorNorm(dir, n, 2);
    if(norm_grad < tg) {
      printf("exit by tg\n");
      printf("%i,\tf(x) %g\t st %g\t|g| %g\n",i, info.function(x1, n), step, norm_grad);
      break;
    }
    step = get_step(step, o_step, stp, info, dir, x, n);
    scaleVector(dir, n, step);
    // printf("alpha*Gradient: ");printVector(dir, n);
    substractVectors(x, dir, x1, n);

    fdif = info.function(x1, n) - fx;
    double dif =0;
    dif = fabs(fdif);
    dif /= fabs(fx) > 1 ? fabs(fx) : 1;

    if( dif < tf ){
      printf("exit by tf\n");
      printf("%i,\tf(x) %g\t st %g\t|g| %g\n",i, info.function(x1, n), step, norm_grad);
      break;
    }
    // printf("%i\t", i);
    double nx = vectorNorm(x, n, 2);

    substractVectors(x, x1, x, n);
    double dif2 = 0;
    dif2 = vectorNorm(x, n, 2);
    dif2 /= nx > 1 ? nx : 1;


    // printf("%i,\tf(x) %g\t st %g\t|g| %g\n",i, info.function(x1, n), step, norm_grad); //printVector(x1, n);
    // fprintf(file, "%lg %lg %lg\n", info.function(x1, n), x1[0], x1[1]);
    // if(i % 1 == 0){
    //   printf("\\multicolumn{1}{|l|}{%i}\t&\t", i);
    //   printf("\\multicolumn{1}{l|}{%lg}\t&\t", dif2);
    //   printf("\\multicolumn{1}{l|}{%lg}\t&\t", vectorNorm(dir, n, 2));
    //   printf("\\multicolumn{1}{l|}{%lg} \\\\ \\hline\n", info.function(x1, n));
    // }
    if( dif2 < tx ){
      printf("exit by tx\n");
      printf("%i,\tf(x) %g\t st %g\t|g| %g\n",i, info.function(x1, n), step, norm_grad);
      break;
    }
    for (int i = 0; i < n; ++i) x[i] = x1[i];
  }
  for (int i = 0; i < n; ++i)
  {
    fprintf(file, "%lg\n", x1[i]);
  }
  // printf("%g\n",info.function(x, n));
  fclose(file);
  free(x1);
  free(dir);
  return info.function(x, n);
}

//rosembrock function evaluator
//n should be at least = 2
double rosembrock(double* x, int n){
  double z = 0.0;
  for(int i = 0; i < n - 1; i++) {
    z += (100.0 * SQUARE((x[i + 1] - SQUARE(x[i]))) + SQUARE(1.0 - x[i]));
  }

  return z;
}
//n should be at least = 2
double* rosembrock_gradient(double* x, int n, double *g){
  g[0] = -400.0 * x[0] * (x[1] - SQUARE(x[0])) - 2.0 * (1.0 - x[0]);
  for (int i = 1; i < (n - 1); ++i)
  {
    g[i] = 200.0 * (x[i] - SQUARE(x[i - 1])) - 400.0 * (x[i+1] - SQUARE(x[i])) * x[i] - 2.0 * (1.0 - x[i]);
  }
  g[n - 1] = 200.0 * (x[n - 1] - SQUARE(x[n - 2]));
  return g;
}
double** rosembrock_hessian(double* x, int n, double **h){
  h[0][0] = -400*(x[1] - 3*SQUARE(x[0])) + 2;
  h[0][1] = -400* x[0];
  for (int i = 2; i < n; ++i) h[0][i] = 0;

  for (int i = 1; i < n - 1; ++i)
  {
    for (int j = 0; j < i - 1; ++j) h[i][j] = 0;

    h[i][i - 1] = -400 * x[i - 1];
    h[i][i] = 202.0 + 1200 * SQUARE(x[i]) - 400.0 * x[i+1];
    h[i][i + 1] = -400* x[i];

    for (int j = i + 2; j < n; ++j) h[i][j] = 0;
  }
  for (int i = 0; i < n-2; ++i) h[n  - 1][i] = 0;
  h[n - 1][n - 2] = -400 * x[n - 2];
  h[n - 1][n - 1] = 200.0;

  return h;
}

// wood's
//keep signature of rosembrock
double wood(double *x, int n){
  double z = 100.0 * SQUARE((SQUARE(x[0]) - x[1]));
  z += SQUARE(x[0] - 1.0) + SQUARE(x[2] - 1.0);
  z += 90.0 * SQUARE(x[2] * x[2] - x[3]);
  z += 10.1 * (SQUARE(x[1] - 1.0) + SQUARE(x[3] - 1.0));
  z += 19.8 * (x[1] - 1.0) * (x[3] - 1.0);
  return z;
}
double* wood_gradient(double *x, int n, double *g){
  g[0] =  400.0 * (SQUARE(x[0]) - x[1]) * x[0] + 2 * (x[0] - 1);
  g[1] = -200.0 * (SQUARE(x[0]) - x[1]) + 20.2 * (x[1]  - 1) + 19.8 * (x[3] - 1);
  g[2] =  2.0   * (x[2] - 1.0)  + 360.0 * (SQUARE(x[2]) - x[3]) * x[2];
  g[3] = -180.0 * (SQUARE(x[2]) - x[3]) + 20.2 * (x[3] - 1.0) + 19.8 * (x[1] - 1.0);
  return g;
}
double** wood_hessian(double *x, int n, double **h){
  // cleanMtx(h, n, n); //zeros mtx
  h[0][0] =  400.0 * (3.0 * SQUARE(x[0]) - x[1]) + 2.0;
  h[2][2] =  2.0 + 360.0 * (3.0 * SQUARE(x[2]) - x[3]);

  h[1][1] =  220.2;
  h[3][3] =  200.2;

  h[0][1] = h[1][0] = -400.0 * x[0];
  h[2][3] = h[3][2] = -360.0 * x[2];
  h[1][3] = h[3][1] = 19.8;
  return h;
}

// smoothing
//to keep signature, y is after x[n] lambda is in x[2n + 1], x should have 2n+1 len
double smoothing(double *x, int n){
  double lambda = x[2 * n];
  double *y = &x[n];
  double z = 0;
  for (int i = 0; i < n - 2; ++i)
  {
    z += SQUARE(x[i] - y[i]) + lambda * SQUARE(x[i + 1] - x[i]);
  }
  z += (x[n - 1] - y[n - 1]);
  return z;
}
double* smoothing_gradient(double *x, int n, double *g){
  double lambda = x[2 * n];
  double *y = &x[n];
  for (int i = 0; i < n; ++i)
  {
    g[i] = 0;
  }

  g[0] = 2 * (x[0] - y[0]) - 2 * lambda * (x[1] - x[0]);

  for (int i = 1; i < n-1; ++i)
  {
    g[i] = 2 * (x[i] - y[i]) + 2 * lambda * (2 * x[i] - x[i - 1] - x[i + 1]);
  }

  g[n-1] = 2 * (x[n-1] - y[n-1]) + 2 * lambda * (x[n-1] - x[n - 2]);
  return g;
}
double** smoothing_hessian(double *x, int n, double **h){
  double lambda = x[2 * n];
  double *y = &x[n];

  h[0][0] = 2 + 2 * lambda;
  h[0][1] = -2 * lambda;

  for (int i = 1; i < n; ++i)
  {
    h[i][i - 1] = -2 * lambda;
    h[i][i]     =  2 + 4 * lambda;
    h[i][i + 1] = -2 * lambda;
  }
  
  h[n - 1][n - 2] = -2 * lambda;
  h[n - 1][n - 1] = 2 + 2 * lambda;

  return h;
}


double cuadraticSolver(double a, double b, double c){
  return (-b + sqrt( SQUARE(b) - 4.0 * a * c)) / (2.0 * a);
}

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

