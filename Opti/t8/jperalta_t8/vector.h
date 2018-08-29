#ifndef VECTOR
#define VECTOR
//memory allocation
double *newVector(int n);
double *newCleanVector(int n);
double *filledVector(int n, double val);

// reset vector
void cleanVector(double* v, int n);
void resetVector(double* v, int n, double val);
//scale vector
void scaleVector(double *v, int n, double val);
void normalizeVector(double *v, int n);
void copyVector(double *v, double *c, int n);
//vector operations
//all have same 'signature' -> input1, intput2, outpit, size
void sumVectors(double *v1, double *v2, double *out, int size);
void substractVectors(double *v1, double *v2, double *out, int size);
void multiplyVectors(double *v1, double *v2, double *out, int size);

double dotproduct(double *v1, double *v2, int n);
//norms
//returns pth norm
double vectorNorm(double *v, int n, int p);
double vectorNormEuclideanSquare(double *v, double n);
double vectorNormInfinity(double *v, double n);
//utils
void printVector(double *v, int n);
#endif