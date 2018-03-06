#ifndef OPTIMIZATION
#define OPTIMIZATION
#define SQUARE(x) ((x)*(x))

typedef enum step { StepFijo, StepAprox, StepHess} Step;

typedef struct fncinfo {
  double (*function)(double*, int n);
  double* (*gradient)(double*, int n, double*);
  double** (*hessian)(double*, int n, double**);
} FuncInfo;


/**
  Set of simple optmization functions that use steepest descent
**/
double get_step_hess(FuncInfo info, double *x, int n, double *g);
double get_step_approx(FuncInfo info, double *x, int n, double alp, double* g);
double optimize_function(FuncInfo info, Step stp, double *x, int n, int iter, double tg, double tx, double tf, double step);


/**
  Set of test functions to optimiza, include function, gradient and hessian
**/

//rosembrock
double rosembrock(double* x, int n);
double* rosembrock_gradient(double* x, int n, double *g);
double** rosembrock_hessian(double* x, int n, double **h);

//wood
double wood(double *x, int n);
double* wood_gradient(double *x, int n, double *g);
double** wood_hessian(double *x, int n, double **h);

//smoothing function
double smoothing(double *x, int n);
double* smoothing_gradient(double *x, int n, double *g);
double** smoothing_hessian(double *x, int n, double **h);
#endif