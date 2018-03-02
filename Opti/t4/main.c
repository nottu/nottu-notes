#include <stdio.h>
#include <float.h>
#include <math.h>
#include "matrix.h"

double Rosenbrock(double *x, int n){
  double aux=0;
  double sum1=0, sum2=0;
  for(int i=0;i<(n-1);i++){
    sum1+=(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]);
    sum2+=(1-x[i]*x[i])*(1-x[i]*x[i]);
  }
  aux=100*sum1 + sum2;
  return(aux);
}

double* gRosenbrock(double *x, int n, double *out){
  for(int i=1;i<(n-1);i++){
    out[i]=-400*x[i]*(x[i+1]-x[i]*x[i])-2*(1-x[i])+200*(x[i]-x[i-1]*x[i-1]);
  }
  out[0]=-400*(x[1]-x[0]*x[0])*x[0]-2*(1-x[0]);
  out[n-1]=200*(x[n-1]-x[n-2]*x[n-2]);
  return out;
}


double** hRosenbrock(double *x, int n, double **out){
  for(int i=1;i<(n-1);i++){
    out[i][i]=202+1200*x[i]*x[i]-400*x[i+1];
    out[i][i-1]=-400*x[i-1];
    out[i-1][i]=-400*x[i-1];
  }
  out[0][0]=1200*x[0]*x[0]-400*x[1]+2;
  out[n-1][n-1]=200; 
  out[n-2][n-1]=-400*x[n-2];
  out[n-1][n-2]=-400*x[n-2];
  return out;
}

typedef enum step { StepFijo, StepAprox, StepHess} Step;

typedef struct fncinfo {
  double (*function)(double*, int n);
  double* (*gradient)(double*, int n, double*);
  double** (*hessian)(double*, int n, double**);
} FuncInfo;

double get_step_hess(FuncInfo info, double *x, int n, double *g){
  double gtg = productoPunto(g, g, n);

  double **h = allocMtx(n, n);
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

double get_step_approx(FuncInfo info, double *x, int n, double alp, double fdif, double* g){
  double gtg = productoPunto(g, g, n);
  alp = gtg * pow(alp, 2) / (2 * (fdif + alp*gtg) );
  return alp;
}

double optimize_function(FuncInfo info, Step stp, double *x, int n, int iter, double tg, double tx, double tf){
  double step = 0;
  double *x1 = (double*)malloc(sizeof(double)*n);
  double *dir = (double*)malloc(sizeof(double)*n);
  double fdif = 0;
  for (int i = 0; i < iter; ++i)
  {
    info.gradient(x, n, dir); //max increse, use negative step...
    if(norma2Vect(dir, n) < tg){
      break;
    }
    switch(stp){
      case StepFijo:
        step = 0.0001;
        break;
      case StepHess:
        step = get_step_hess(info, x, n, dir);
      case StepAprox:
        if(i == 0) step = 0.0005;
        else step = get_step_approx(info, x, n, step, fdif, dir);
        break;
      default:
        break;
    }
    // normalizaVect(dir, n);
    //printVect(x, 2);
    vectorScalar(dir, step, n);
    restaVector(x, dir, x1, n);
    printf("%i, step %g, %g\n",i, step, info.function(x1, n));// printVect(x1, n);
    double fx = info.function(x, n);
    fdif = info.function(x1, n) - fx;
    double dif = fabs(fdif);
    dif /= fx > 1 ? fx : 1;
    if( dif < tf ){
      break;
    }
    double nx = norma2Vect(x, n);
    restaVector(x, x1, x, n);
    dif = norma2Vect(x, n);
    dif /= nx > 1 ? nx : 1;

    if( dif < tf ){
      break;
    }

    double *xt = x;
    x = x1;
    x1 = xt;
  }
  free(dir);
  return info.function(x1, n);
}

//rosembrock function evaluator
//n should be at least = 2
double rosembrock(double* x, int n){
  double z = 0;
  for (int i = 0; i < n - 1; ++i)
  {
    z += 100 * pow((x[i+1] - pow(x[i], 2)) , 2)  +  pow((1 - x[i]), 2);
  }
  return z;
}
//n should be at least = 2
double* rosembrock_gradient(double* x, int n, double *g){
  g[0] = -400 * (x[1] - pow(x[0], 2)) * x[0] - 2 * (1 - x[0]);
  for (int i = 1; i < (n - 1); ++i)
  {
    g[i] =200*(x[i] - pow(x[i-1], 2)) - 400 * (x[i+1] - pow(x[i], 2)) * x[i] - 2 * (1 - x[i]);
  }
  g[n - 1] = 200 * (x[n-1] - pow(x[n-2], 2));
  return g;
}

double** rosembrock_hessian(double* x, int n, double **h){
  h[0][0] = -400*(x[1] - 3*pow(x[0], 2)) + 2;
  h[0][1] = -400* x[0];
  for (int i = 2; i < n; ++i) h[0][i] = 0;

  for (int i = 1; i < n - 1; ++i)
  {
    for (int j = 0; j < i - 1; ++j) h[i][j] = 0;

    h[i][i - 1] = -400* x[i - 1];
    h[i][i] = 202 + 1200 * pow(x[i], 2) - 400 * x[i+1];
    h[i][i + 1] = -400* x[i];

    for (int j = i + 2; j < n; ++j) h[i][j] = 0;
  }
  for (int i = 0; i < n-2; ++i) h[n  - 1][i] = 0;
  h[n - 1][n - 2] = -400 * x[n - 2];
  h[n - 1][n - 1] = 200;

  return h;
}

// wood's
//keep signature of rosembrock
double wood(double *x, int n){
  double z = 100 * pow(pow(x[0], 2) - x[1], 2);
  z += pow(x[0] - 1, 2) + pow(x[2] - 1, 2);
  z += 90 * pow(pow(x[2], 2) - x[3], 2);
  z += 10.1 * (pow(x[1] - 1, 2) + pow(x[3] - 1, 2));
  z += 19.8 * (x[1] - 1) * (x[3] - 1);
  return z;
}

double *wood_gradient(double *x, int n, double *g){
  g[0] = 400 * (pow(x[0], 2) - x[1]) * x[0] + 2 * (x[0] - 1);
  g[1] = -200 * (pow(x[0], 2) - x[1]) + 20.2 * (x[1] - 1) + 19.8 * (x[3] - 1);
  g[2] = 2 * (x[2] - 1) + 360 * (pow(x[2], 2) - x[3]) * x[2];
  g[3] = -180 * (pow(x[2], 2) - x[3]) + 20.2 * (x[3] - 1) + 19.8 * (x[1 - 1]);

  return g;
}

double ** wood_hessian(double *x, int n, double **h){
  cleanMtx(h, n, n); //zeros mtx
  h[0][0] = 400 * (3 * pow(x[0], 2) - x[1]) + 2;
  h[0][1] = -400 * x[0];
  h[1][0] = -400 * x[0];

  h[1][1] = 220.2;
  h[1][3] = -19.8;
  h[2][2] = 2 + 360 * (3 * pow(x[2], 2) - x[3]);
  h[2][3] = -360*x[2];
  h[3][3] = 200.2;
  h[3][2] = -360*x[2];
  h[3][1] = -19.8;
  return h;
}

int main(int argc, char const *argv[]){

  FuncInfo f;
  int n;
  n = 2;
  f.function = Rosenbrock;
  f.gradient = gRosenbrock;
  f.hessian  = hRosenbrock;
  // f.function = rosembrock;
  // f.gradient = rosembrock_gradient;
  // f.hessian  = rosembrock_hessian;
  double x[] = {-1.2, 1};
  // n = 100;
  // f.function = rosembrock;
  // f.gradient = rosembrock_gradient;
  // f.hessian  = rosembrock_hessian;
  // double x[n];
  // for (int i = 0; i < n; ++i) x[i] = 1;
  // x[0] = -1.2;
  // x[98] = -1.2;
  // n = 4;
  // f.function = wood;
  // f.gradient = wood_gradient;
  // f.hessian = wood_hessian;
  // double x[] = {-3, -1, -3, -1};

  optimize_function(f, StepHess, x, n, 100000, 1e-6, 1e-12, 1e-12);

  return 1;
}