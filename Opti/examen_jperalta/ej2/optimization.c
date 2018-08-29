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
    scaleVector(g, n, alpha);         // alpha d
    sumVectors(x, g, g, n);           //x + alpha d
    double fxk = info.function(g, n); //f(xk)
    if(fxk <= fx + alpha*c*fd) {
      // printf("fxk %g\n", fxk);
      break;
    }
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
  double o_step = step;
  for (int i = 0; i < iter; ++i)
  {
    fx = info.function(x, n);
    info.gradient(x, n, dir); //max increse, use negative step...
    double norm_grad = vectorNorm(dir, n, 2);
    if(norm_grad < tg) {
      // printf("exit by tg\n");
      // printf("%i,\tf(x) %g\t st %g\t|g| %g\n",i, info.function(x1, n), step, norm_grad);
      break;
    }
    step = get_step(step, o_step, stp, info, dir, x, n);
    // printf("Dir := ");printVector(dir, n);
    // printf("stp = %g\n", step);
    scaleVector(dir, n, step);
    substractVectors(x, dir, x1, n);
    fdif = info.function(x1, n) - fx;
    double dif =0;
    dif = fabs(fdif);
    dif /= fabs(fx) > 1 ? fabs(fx) : 1;

    // if( dif < tf ){
    //   printf("exit by tf\n");
    //   printf("%i,\tf(x) %g\t st %g\t|g| %g\n",i, info.function(x1, n), step, norm_grad);
    //   break;
    // }

    double nx = vectorNorm(x, n, 2);

    substractVectors(x, x1, x, n);
    double dif2 = 0;
    dif2 = vectorNorm(x, n, 2);
    dif2 /= nx > 1 ? nx : 1;

    // if( dif2 < tx ){
    //   printf("exit by tx\n");
    //   printf("%i,\tf(x) %g\t st %g\t|g| %g\n",i, info.function(x1, n), step, norm_grad);
    //   break;
    // }
    for (int i = 0; i < n; ++i) x[i] = x1[i];
    // printVector(x, n);
  }
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
double cauchy_point(int n, double *g, double **h){
  double gtg = dotproduct(g, g, n);

  double *hg = newVector(n);

  multiplyMatrixVector(h, g, n, n, hg);
  double alp = gtg / dotproduct(g, hg, n);

  free(hg);
  return -1 * alp;
}

void doglegDirection(int n, double *g, double **h, double* cauchy_dir, double reg_sz, double *dir){
  double *pb = luSolver2(h, g, n, n);
  if(pb == NULL) return;
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

double taylorEval(double fx, int n, double *g, double **h, double *p){
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

  double fx;
  double norm2 = vectorNorm(g, n, 2);
  double reg_sz = reg_szM;
  int iter = 0;

  double *pu = newVector(n);
  double *dir = newVector(n);
  double *x1 = newVector(n);
  while(norm2 >= tg && iter < max_iter){
    fx = info.function(x, n);
    iter++;
    double a = cauchy_point(n, g, H);
    copyVector(g, pu, n);
    scaleVector(pu, n, a);
    double pu_norm = vectorNorm(pu, n ,2);
    copyVector(pu, dir, n);
    if(pu_norm < reg_sz){
      if(isdefpos(H, n)) {
        //do dogleg direction
        // printf("DOGLEG!\n");
        doglegDirection(n, g, H, pu, reg_sz, dir);
      }
    } else {
      //make smaller
      scaleVector(dir, n, reg_sz/pu_norm);
    }

    sumVectors(x, dir, x1, n); //x = x + dir
    //check if solution is good!
    double taylorev = fx - taylorEval(fx, n, g, H, dir);
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
  free(pu);
  free(dir);
  freeMatrix(H);
  return fx;
}

void linealConjugateGradient(FuncInfo info, double *x, int n, double tg, double tx, double tf){
  double *g  = info.gradient(x, n, newVector(n));
  double **h = info.hessian(x, n, allocMatrix(n, n));
  double *d = newVector(n);
  copyVector(g, d, n);
  scaleVector(d, n, -1);

  double *hg = newVector(n);
  int iter = 0;
  printf("f(x): %g\t||g(x)||: %g\n", info.function(x, n), vectorNorm(g, n ,2));
  // printf("Gradient\t");printVector(g, n);
  while(vectorNorm(g, n ,2) > tg && iter++ < n){
    multiplyMatrixVector(h, d, n, n, hg);
    double dQd = dotproduct(d, hg, n);
    double alpha = -1 * dotproduct(g, d, n) / dQd;

    copyVector(d, g, n);
    scaleVector(g, n, alpha); //alpha_k dk
    if(vectorNorm(g, n ,2) < tx) break;

    // printVector(x, n);
    sumVectors(x, g, x, n); // x = x + alpha_k d_k
    // printVector(x, n);

    g = info.gradient(x, n, g);
    h = info.hessian(x, n, h);

    multiplyMatrixVector(h, g, n, n, hg);
    double beta = dotproduct(d, hg, n) / dQd;
    scaleVector(d, n, beta);
    substractVectors(d, g, d, n);
    printf("Iter: %i, f(x): %g\t||g(x)||: %g\t alpha %g\n", iter, info.function(x, n), vectorNorm(g, n ,2), alpha);
    if(info.function(x, n) < tf) break;
  }

  printf("X*\t");printVector(x, n);

  free(d);
  free(g);
  freeMatrix(h);
}
double fletcher_reeves(double *g, double *g1, double *d, int n){
  return dotproduct(g1, g1, n) / dotproduct(g, g, n);
}
double polak_ribiere(double *g, double *g1, double *d, int n){
  double *dif = newVector(n);
  substractVectors(g1, g, dif, n);
  double beta = dotproduct(g1, dif, n) / dotproduct(g, g, n);
  free(dif);
  return beta > 0 ? beta : 0; //max(0, beta)
}
double hestenes_stiefel(double *g, double *g1, double *d, int n){
  double *dif = newVector(n);
  substractVectors(g1, g, dif, n);
  double beta = dotproduct(g1, dif, n) / dotproduct(dif, d, n);
  free(dif);
  return beta;
}

void non_linear_conjugateGradient(FuncInfo info, double *x, int n, double tg, double tx, double tf, BetaFun betaFun){
  double *g  = info.gradient(x, n, newVector(n));
  double *d = newVector(n);
  copyVector(g, d, n);
  scaleVector(d, n, -1);

  double *ad = newVector(n);
  double *g1 = newVector(n);
  int iter = 0;
  FILE *out = fopen("out.txt", "w");
  printf("f(x): %g\t||g(x)||: %g\n", info.function(x, n), vectorNorm(g, n ,2));
  fprintf(out, "%g \t %g \t %g\n", info.function(x, n), vectorNorm(d, n ,2), vectorNorm(g, n ,2));
  while(vectorNorm(g, n ,2) > tg && iter++ < 100){
    double alpha = get_step_backtracking(info, d, x, 1E-4, 0.95, 1, n); //could use other method...

    copyVector(d, ad, n);
    scaleVector(ad, n, alpha); //alpha_k dk
    sumVectors(x, ad, x, n); // x = x + alpha_k d_k
    // if(alpha < 1E-10) break;
    if(vectorNorm(ad, n ,2) < tx) break;

    g1 = info.gradient(x, n, g1);
    double beta = betaFun(g, g1, d, n);
    scaleVector(d, n, beta);
    copyVector(g1, g, n);
    substractVectors(d, g, d, n);
    printf("Iter: %i, f(x): %g\t||g(x)||: %g\t alpha %g\t beta %g\n", iter, info.function(x, n), vectorNorm(g, n ,2), alpha, beta);
    fprintf(out, "%g \t %g \t %g\n", info.function(x, n), vectorNorm(d, n ,2), vectorNorm(g, n ,2));
    if(info.function(x, n) < tf) break;
  }
  fclose(out);
  free(g);
  free(d);
  free(g1);
  free(ad);
}

//convex 1
double convex1(double *x, int n){
  double val = 0;
  for (int i = 0; i < n; ++i) {
    val += exp(x[i]) - x[i];
  }
  return val;
}
double* convex1_gradient(double *x, int n, double *g){
  for (int i = 0; i < n; ++i) {
    g[i] = exp(x[i]) - x[i];
  }
  return g;
}
double** convex1_hessian(double *x, int n, double **h){
  return NULL;
}

//convex 2
double convex2(double *x, int n){
  double val = 0;
  for (int i = 0; i < n; ++i) {
    val += i * (exp(x[i]) - x[i]) / n;
  }
  return val;
}
double* convex2_gradient(double *x, int n, double *g){
  for (int i = 0; i < n; ++i) {
    g[i] = i * (exp(x[i]) - x[i]) / n;
  }
  return g;
}
double** convex2_hessian(double *x, int n, double **h){
  return NULL;
}

void printInfo_nolineal(int k, double *x, double norm, double cond, int n){
  printf("%d\t", k);
  printf("\t||F|| = %.7g\tK(J) = %g\tx = (", norm, cond);
  printVector(x, n);
  printf("\b)\n");
}

void printInfo_fun(int k, double *x, double norm, double cond, int n){
  printf("%d\t", k);
  printf("\tf(x) = %.7g\tK(H) = %g\tx = (", norm, cond);
  printVector(x, n);
  printf("\b)\n");
}

double *newton_no_lineal(double *x, MultiFun f, Jacobian j, int max_iter, double toler, int n){
  double *fx  = f(x, n, newVector(n));
  double *x1 = newVector(n);
  double **jx = j(x, n, allocMatrix(n, n));

  copyVector(x, x1, n);
  int iter = 0;
  printInfo_nolineal(iter, x1, vectorNorm(fx, n, 2), matrix_condition(jx, n), n);
  while(vectorNorm(fx, n, 2) > toler && iter++ < max_iter){
    scaleVector(fx, n, -1);
    double *x_d = luSolver2(jx, fx, n, n);

    for (int i = 0; i < n; ++i){
      x1[i] += x_d[i];
    }
    free(x_d);
    f(x1, n, fx);
    j(x1, n, jx);
    printInfo_nolineal(iter, x1, vectorNorm(fx, n, 2), matrix_condition(jx, n), n);
  }
  free(fx);
  freeMatrix(jx);
  return x1;
}

double *broyden(double *x, MultiFun f, Jacobian j, int max_iter, double toler, int n){
  double *fx  = f(x, n, newVector(n));
  double *x1 = newVector(n);
  double *f_delta = newVector(n);
  double *f_aux = newVector(n);

  double **jx = j(x, n, allocMatrix(n, n));
  double **jx_d = allocMatrix(n, n);

  copyVector(x, x1, n);
  int iter = 0;
  printInfo_nolineal(iter, x1, vectorNorm(fx, n, 2), matrix_condition(jx, n), n);
  while(vectorNorm(fx, n, 2) > toler && iter++ < max_iter){
    copyVector(fx, f_delta, n);
    scaleVector(fx, n, -1);
    double *x_d = luSolver2(jx, fx, n, n);

    for (int i = 0; i < n; ++i){
      x1[i] += x_d[i];
    }
    f(x1, n, fx);

    printInfo_nolineal(iter, x1, vectorNorm(fx, n, 2), matrix_condition(jx, n), n);

    substractVectors(fx, f_delta, f_delta, n);
    double x_d_n = dotproduct(x_d, x_d, n);
    multiplyMatrixVector(jx, x_d, n, n, f_aux);
    substractVectors(f_delta, f_aux, f_delta, n);
    scaleVector(f_delta, n, 1 / x_d_n);
    for (int i = 0; i < n; ++i){
      for (int j = 0; j < n; ++j){
        jx_d[i][j] = f_delta[i] * x_d[j];
      }
    }
    addMatrix(jx, jx_d, n, n, jx);
    free(x_d);
  }
  free(fx);
  free(f_delta);
  freeMatrix(jx);
  freeMatrix(jx_d);
  return x1;
}
double **HessianAproximator(double *x, FuncInfo info, int n, double h){
  double **A = allocMatrix(n, n);
  double **I = allocMatrixIdentity(n, n);
  scaleMatrix(I, n, n, h);
  double *e = newVector(n);
  double fx = info.function(x, n);
  double hh = SQUARE(h);
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      A[i][j] = fx;
      sumVectors(x, I[i], e, n);
      A[i][j] -= info.function(e, n);
      sumVectors(e, I[j], e, n);
      A[i][j] += info.function(e, n);
      sumVectors(x, I[j], e, n);
      A[i][j] -= info.function(e, n);
      A[i][j] /= hh;
      A[j][i] = A[i][j];
    }
  }
  free(e);
  freeMatrix(I);
  return A;
}

double *bfgs(double *x, FuncInfo info, double **H, int max_iter, double toler, int n){
  double *x1 = newVector(n);
  copyVector(x, x1, n);
  double fx  = info.function(x1, n);
  double *g  = info.gradient(x1, n, newVector(n));
  double *g2 = newVector(n);
  double *p  = newVector(n);
  double *s  = newVector(n);

  double **rsy = allocMatrix(n, n);
  double **rss = allocMatrix(n, n);
  double **I = allocMatrixIdentity(n, n);
  double **H2 = allocMatrix(n, n);
  int iter = 0;
  printInfo_fun(iter, x1, fx, matrix_condition(H, n), n);
  while(vectorNorm(g, n, 2) > toler && iter < max_iter ){
    multiplyMatrixVector(H, g, n, n, p);
    scaleVector(p, n, -1);

    copyVector(p, s, n);
    double a = get_step_backtracking(info, p, x1, 1E-4, 0.9, 1.0, n);
    scaleVector(s, n, a);
    sumVectors(x1, s, x1, n);

    copyVector(g, g2, n);
    g = info.gradient(x1, n, g);
    substractVectors(g, g2, g2, n); //g2 => y
    double r = dotproduct(g2, s, n);
    // printf("r %g\n", r);
    if(r > 1E-4){
      r = 1 / r;
      vectorVector(s, g2, n, n, rsy);//sy
      scaleMatrix(rsy, n, n, r);//rsy
      substractMatrix(I, rsy, n, n, rsy);//I - rsy

      multiplyMatrix(rsy, H, n, n, n, H2);

      vectorVector(g2, s, n, n, rsy);//ys
      scaleMatrix(rsy, n, n, r);//rsy
      substractMatrix(I, rsy, n, n, rsy);//I - rys

      multiplyMatrix(H2, rsy, n, n, n, H);

      vectorVector(s, s, n, n, rss);
      scaleMatrix(rss, n, n, r);

      addMatrix(H, rss, n, n, H);
      // makedefpos(H, n);
    }
    double dif = fx;
    fx  = info.function(x1, n);
    dif -= fx;
    iter++;
    // printMatrix(H, n, n);printf("\n");
    printInfo_fun(iter, x1, fx, matrix_condition(H, n), n);
    if(dif < 1E-6) break;
  }

  free(g);
  free(g2);
  free(p);
  free(s);
  freeMatrix(rsy);
  freeMatrix(rss);
  freeMatrix(I);
  freeMatrix(H2);

  return x1;
}

double *levenbergMarquardt(double *x_0, MultiFun f, Jacobian2 j, int max_iter, double toler, int n, int m, double v){
  double *x1 = newVector(n);
  double *xn = newVector(n);
  copyVector(x_0, x1, n);
  double *r  = f(x1, m, newVector(m));
  double **J = j(x1, n, m, allocMatrix(m, n));
  double **JJ = allocMatrix(n, n);
  multiplyMatrixTransposed(J, J, n, n, m, JJ);
  double l = 0;
  for (int i = 0; i < n; ++i) l = MAX(l, JJ[i][i]);
  double **lI = allocMatrixIdentity(n, n);
  double *g = multiplyMatrixTransposedVector(J, r, m, n, newVector(n));
  int i = 0;
  double norm = vectorNorm(g, n, 2);
  printf("%i\t%g, %g\t%g\t%g\t%g\n", i, x1[0], x1[1], l, dotproduct(r, r, m)/2, norm);
  while(norm > toler && i < max_iter){
    for (int j = 0; j < n; ++j) lI[j][j] = l;
    addMatrix(JJ, lI, n, n, JJ);
    scaleVector(g, n, -1);
    double *p = luSolver2(JJ, g, n, n);
    double fx = dotproduct(r, r, m); //should divide by 2
    sumVectors(x1, p, xn, n);// xn = x1 + p
    free(p);
    f(xn, m, r); //r
    double fxn = dotproduct(r, r, m); //should divide by 2
    if(fxn <= fx){
      copyVector(xn, x1, n);
      l /= v;
      j(x1, n, m, J); //recalculate J
      multiplyMatrixTransposed(J, J, n, n, m, JJ);// recalculate JJ
      multiplyMatrixTransposedVector(J, r, m, n, g);//recalc g
      norm = vectorNorm(g, n, 2);
    } else {
      l *= v;
      r = f(x1, m, r); //reset r
      scaleVector(g, n, -1); //reset g
    }
    // printf("%i\t%g, %g\t%g\t%g\t%g\n", i, x1[0], x1[1], l, dotproduct(r, r, m)/2, norm);
    i++;
  }
  printf("%i\t%g, %g\t%g\t%g\t%g\n", i, x1[0], x1[1], l, dotproduct(r, r, m)/2, norm);
  printf("\n");
  free(g);
  free(r);
  free(xn);
  freeMatrix(lI);
  freeMatrix(JJ);
  freeMatrix(J);
  return x1;
}
double* finete_difference(double *x, FunEval f, int n, double h, double *g){
  double *x_d = newVector(n);
  copyVector(x, x_d, n);
  double eval = f(x, n);
  for (int i = 0; i < n; ++i){
    x_d[i] += h;
    double e2 = f(x_d, n);
    g[i] = e2 - eval;
    x_d[i] -= h;
  }
  free(x_d);
  return g;
}
double* quadraticPenalizaition(double mu, double *x, double toler, int max_iter, int n,
                                FunEval f, FunEval *g, FunEval *h, int m, int l){
  // en gcc podemos definir funciones dentro de funciones...
  // printVector(x, n);printf("\n");
  double q_x(double *x, int n){
    double res = f(x, n);
    for (int i = 0; i < m; ++i){
      double v = g[i](x, n);
      res += (mu/2.0) * SQUARE(v);
    }
    for (int i = 0; i < l; ++i){
      double v = h[i](x, n);
      v = MAX(0, v);
      res += (mu/2.0) * SQUARE(v);
    }
    return res;
  }
  double *q_x_g(double *x, int n, double *grad){return finete_difference(x, q_x, n, 0.01, grad);}
  FuncInfo info; info.function = q_x; info.gradient = q_x_g;

  double *x_1 = newVector(n); copyVector(x, x_1, n);
  int i = 0;
  for (i = 0; i < max_iter; ++i){
    // printVector(x_1, n); printf("\n");
    double t = 0.1/(i+1);
    double v = optimize_function(info, StepBacktrack, x_1, n, 1000, toler, t, toler, 1.0);
    double p = 2.0 * (q_x(x_1, n) - f(x_1, n)) / mu;
    // printf("p = %g\n", p);
    if( p < toler ) break;
    mu *= 3.0;
  }
  for (int j = 0; j < m; ++j){
    printf("m_%i : %g\t", j+1, g[j](x_1, n));
  }
  for (int j = 0; j < l; ++j){
    printf("h_%i : %g\t", j+1, h[j](x_1, n));
  }
  printf("\n");
  return x_1;
}