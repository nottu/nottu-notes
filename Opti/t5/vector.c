#include "vector.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//memory allocation
double *newVector(int n){
  double *v = (double*)malloc(sizeof(double) * n);
  return v;
}
double *newCleanVector(int n){
  return filledVector(n, 0);
}
double *filledVector(int n, double val){
  double *v = newVector(n);
  resetVector(v, n, val);
  return v;
}
// reset vector
void cleanVector(double* v, int n){
  resetVector(v, n, 0);
}
void resetVector(double* v, int n, double val){
  for (int i = 0; i < n; ++i) v[i] = val;
}
//scale vector
void scaleVector(double *v, int n, double val){
  for (int i = 0; i < n; ++i) v[i] *= val;
}
//euclidean normalization
void normalizeVector(double *v, int n){
  double norm = vectorNorm(v, n, 2); //find norm
  scaleVector(v, n, 1.0/norm);
}
//vector operations
//all have same 'signature' -> input1, intput2, outpit, size
void sumVectors(double *v1, double *v2, double *out, int n){
  for (int i = 0; i < n; ++i) out[i] = v1[i] + v2[i];
}
void substractVectors(double *v1, double *v2, double *out, int n){
  for (int i = 0; i < n; ++i) out[i] = v1[i] - v2[i];
}
void multiplyVectors(double *v1, double *v2, double *out, int n){
  for (int i = 0; i < n; ++i) out[i] = v1[i] * v2[i];
}
double dotproduct(double *v1, double *v2, int n){
  double c = 0;
  for (int i = 0; i < n; i++) {
    c += v1[i] * v2[i];
  }
  return c;
}
//norms
double pthPower(double val, int p){
  double o = 1;
  for (int i = 0; i < p; ++i)
  {
    o *= val;
  }
  return o;
}
//returns pth norm
double vectorNorm(double *v, int n, int p){
  double norm = 0;
  for (int i = 0; i < n; ++i)
  {
    norm+= pthPower(fabs(v[i]), p);
  }
  if( p > 1 ) norm = pow(norm, 1.0/p);
  return norm;
}
double vectorNormEuclideanSquare(double *v, double n){
  double norm = 0;
  for (int i = 0; i < n; ++i)
  {
    norm+= pthPower(fabs(v[i]), 2);
  }
  return norm;
}
double vectorNormInfinity(double *v, double n){
  double max = 0;
  for (int i = 0; i < n; ++i)
  {
    double vabs = fabs(v[i]);
    if(vabs > max) max = vabs;
  }
  return max;
}

//utils
void printVector(double *v, int n){
  for (int i = 0; i < n; ++i) printf("%g ", v[i]);
  printf("\n");
}