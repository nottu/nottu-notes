#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "matrix.h"
/**
Box Muller transform
https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
**/
double generateGaussianNoise(double mu, double sigma)
{
  static const double two_pi = 2.0*M_PI;

  double u1, u2;
  do
   {
     u1 = rand() * (1.0 / RAND_MAX);
     u2 = rand() * (1.0 / RAND_MAX);
   }
  while ( u1 <= DBL_EPSILON );

  double z0;
  z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
  return z0 * sigma + mu;
}

double getGx(double x){
  double g = sin(x)/x + generateGaussianNoise(0, 0.1);
  return g;
}

double* getParts(double min, double max, int n){
  double size = (max - min) / (n-1);
  double *parts = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; ++i)
  {
    parts[i] = min + i * size;
  }
  return parts;
}

//returns n*m matrix for least squares
double** getLSMatrix(int n, double m, double* x){
  double **mtx = allocMtx(n, n);
  double *xij = malloc(sizeof(double) * n);
  // xij[0] = n;
  // for (int i = 1; i < n; ++i)
  // {
  //   double val = 0;
  //   for (int j = 0; j < m; ++j)
  //   {
  //     val += pow(x[j], i);
  //   }
  //   xij[i] = val;
  // }
  for (int i = 0; i < n; ++i)
  {
    // mtx[i][0] = xij[i];
    for (int j = 0; j < n; ++j)
    {
      mtx[i][j] = 0;
      for (int k = 0; k < m; ++k)
      {
        mtx[i][j] += pow(x[k], j+i);
      }
    }
  }
  return mtx;
}
double* getLSVec(int n, int m, double *g, double *x){
  double *vec = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; ++i)
  {
    vec[i] = 0;
    for (int j = 0; j < m; ++j)
    {
      vec[i] += pow(x[j], i) * g[j];
    }
  }
  return vec;
}

double evalPoli(double x, double* a, int n){
  double y = a[0];
  for (int i = 1; i < n; ++i)
  {
    y += a[i]*pow(x, i);
  }
  return y;
}

int main(int argc, char** argv){
  int nparts = atoi(argv[1]);
  int nPol = atoi(argv[2]);
  double *x = getParts(0.1, 10, nparts);
  double *g = (double*)malloc(sizeof(double) * nparts);
  for (int i = 0; i < nparts; ++i) g[i] = getGx(x[i]);
  // printVect(x, nparts);
  // printf("-------------\n");
  // printVect(g, nparts);
  // printf("-------------\n");
  double **mtx = getLSMatrix(nPol, nparts, x);
  double *vec = getLSVec(nPol, nparts, g, x);
  // printMtx(mtx, nPol, nPol);
  // printf("-------------\n");
  // printVect(vec, nPol);
  // printf("-------------\n");
  luFactor2(mtx, nPol, nPol);
  double *a = luSolver2(mtx, vec, nPol, nPol);
  printVect(a, nPol);
  char name[30];
  sprintf(name, "out_%i.txt", nPol);
  FILE *f = fopen(name, "w");
  for (int i = 0; i < nparts; ++i)
  {
    fprintf(f, "%lf %lf %lf\n", x[i], g[i], evalPoli(x[i], a, nPol));
  }
  fclose(f);
  free(vec);
  freeMtx(mtx);
  free(g);
  free(x);
  return 0;
}