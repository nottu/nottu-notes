#include "optimization.h"
#include "matrix.h"

double get_step_hess(FuncInfo info, double *x, int n, double *g){
  double gtg = productoPunto(g, g, n);

  double **h = allocMtxClean(n, n);
  info.hessian(x, n, h);
  
  double *hg = (double*)malloc(sizeof(double) * n);
  multMatrizVect(h, g, n, n, hg);
  
  double alp = gtg / productoPunto(g, hg, n);
  
  if(!isdefpos(h, n)) alp *= -1;
  // makedefpos(h, n);

  free(hg);
  freeMtx(h);
  return alp;
}

double get_step_approx(FuncInfo info, double *x, int n, double alp, double* g){
  double gtg = productoPunto(g, g, n);
  double numerator = gtg * alp * alp;

  double *dir = (double*)malloc(sizeof(double) * n);
  info.gradient(x, n, dir);
  double *x1  = (double*)malloc(sizeof(double) * n);
  vectorScalar(dir, alp, n);
  restaVector(x, dir, x1, n);

  double fx  = info.function(x1, n);
  double pfx = info.function(x,  n);

  double denominator = 2.0 * (fx - pfx + alp * gtg);
  free(dir);
  free(x1);
  return numerator/denominator;
  // return 0.005;
}

double optimize_function(FuncInfo info, Step stp, double *x, int n, int iter, double tg, double tx, double tf, double step){
  double *x1  = (double*)malloc(sizeof(double) * n);
  double *dir = (double*)malloc(sizeof(double) * n);
  double fdif = 0;
  double fx = 0;
  FILE *file = fopen("data.txt", "w");
  // fprintf(file, "%lg %lg %lg\n", info.function(x, n), x[0], x[1]);
  for (int i = 0; i < iter; ++i)
  {
    fx = info.function(x, n);
    info.gradient(x, n, dir); //max increse, use negative step...
    //printVect(dir, n);
    if(norma2Vect(dir, n) < tg) {
      printf("exit by tg\n");
      break;
    }
    switch(stp){
      case StepFijo:
        break;
      case StepHess:
        step = get_step_hess(info, x, n, dir);
        break;
      case StepAprox:
        step = get_step_approx(info, x, n, step, dir);
        break;
      default:
        break;
    }
    vectorScalar(dir, step, n);
    restaVector(x, dir, x1, n);

    fdif = info.function(x1, n) - fx;
    double dif =0;
    dif = fabs(fdif);
    dif /= fabs(fx) > 1 ? fabs(fx) : 1;

    if( dif < tf ){
      printf("exit by tf\n");
      break;
    }
    // printf("%i\t", i);
    double nx = norma2Vect(x, n);

    restaVector(x, x1, x, n);
    double dif2 = 0;
    dif2 = norma2Vect(x, n);
    dif2 /= nx > 1 ? nx : 1;


    printf("%i,\t%g\n",i, info.function(x1, n));// printVect(x1, n);
    // fprintf(file, "%lg %lg %lg\n", info.function(x1, n), x1[0], x1[1]);
    // if(i % 1 == 0){
    //   printf("\\multicolumn{1}{|l|}{%i}\t&\t", i);
    //   printf("\\multicolumn{1}{l|}{%lg}\t&\t", dif2);
    //   printf("\\multicolumn{1}{l|}{%lg}\t&\t", norma2Vect(dir, n));
    //   printf("\\multicolumn{1}{l|}{%lg} \\\\ \\hline\n", info.function(x1, n));
    // }

    if( dif2 < tx ){
      printf("exit by tx\n");
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