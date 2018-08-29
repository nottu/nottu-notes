#ifndef OPTIMIZATION
#define OPTIMIZATION
#define SQUARE(x) ((x)*(x))

#include "matrix.h"
#include "matrix_factor.h"

typedef enum step { StepFijo, StepAprox, StepHess, StepBacktrack, StepInterpol} Step;
typedef double (*BetaFun)(double*, double*, double*, int);

typedef double* (MultiFun)(double *, int, double*);
typedef double** (Jacobian)(double *, int, double**);

typedef struct fncinfo {
  double (*function)(double*, int n);
  double* (*gradient)(double*, int n, double*);
  double** (*hessian)(double*, int n, double**);
} FuncInfo;


/**
  Set of simple optmization functions that use steepest descent
**/
double get_step(double alp, double alp_0, Step stp, FuncInfo info, double *gradient, double *x, int n);
double get_step_hess(FuncInfo info, double *x, int n, double *g);
double get_step_approx(FuncInfo info, double *x, int n, double alp, double* g);
double get_step_backtracking(FuncInfo info, double* d, double *x, double c, double p, double alpha, int n);
double get_step_interpolated(FuncInfo info, double* d, double *x, double c, double alpha, int n);
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

//

double digitEstimator(double *x , int n);
double* digitEstimator_gradient(double *x , int n, double *g);
double** digitEstimator_hessian(double *x , int n, double **h);

//convex 1
double convex1(double *x, int n);
double* convex1_gradient(double *x, int n, double *g);
double** convex1_hessian(double *x, int n, double **h);

//convex 2
double convex2(double *x, int n);
double* convex2_gradient(double *x, int n, double *g);
double** convex2_hessian(double *x, int n, double **h);

//
double cuadraticSolver(double a, double b, double c);

//dogleg
double cauchy_point(int n, double *g, double **h);
void   doglegDirection(int n, double *g, double **h, double* cauchy_dir, double reg_sz, double *dir);
double taylorEval(double fx, int n, double *g, double **h, double *p);
double doglegOptimize(FuncInfo info, double *x, int n, int max_iter, double tg, double tx, double tf, double reg_szM);

//
void linealConjugateGradient(FuncInfo info, double *x, int n, double tg, double tx, double tf);
double fletcher_reeves(double *g, double *g1, double *d, int n);
double polak_ribiere(double *g, double *g1, double *d, int n);
double hestenes_stiefel(double *g, double *g1, double *d, int n);
void non_linear_conjugateGradient(FuncInfo info, double *x, int n, double tg, double tx, double tf, BetaFun betaFun);

//
double *newton_no_lineal(double *x, MultiFun f, Jacobian j, int max_iter, double toler, int n);
double *broyden(double *x, MultiFun f, Jacobian j, int max_iter, double toler, int n);
#endif