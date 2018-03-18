#include "matrix.h"
#include "vector.h"
#include <stdlib.h>
//memory management
double** allocMatrix(int nr, int nc){
  double **mtx = malloc((sizeof(double*) * nr) + sizeof(int));
  int *indi = (int*) mtx;
  mtx = (void*) indi + sizeof(int);
  if(nr * nc * sizeof(double) < MTXMAXSIZE) {
    indi[0] = 0; //indicate 1 block
    mtx[0] = (double*)malloc(sizeof(double) * nr * nc);
    for (int i = 1; i < nr; ++i) {
      mtx[i] = mtx[i-1] + nc;
    }
  } else {
    indi[0] = nr; //indicate nr block
    for (int i = 0; i < nr; ++i) {
      mtx[i] = malloc(sizeof(double) * nc);
    }
  }
  return mtx;
}
double** allocMatrixClean(int nr, int nc){
  double **mtx = allocMatrix(nr, nc);
  for (int i = 0; i < nr; ++i)
    for (int j = 0; j < nc; ++j)
      mtx[i][j] = 0;
  return mtx;
}
double** allocMatrixIdentity(int nr, int nc){
  double **mtx = allocMatrix(nr, nc);
  for (int i = 0; i < nr; ++i){
    for (int j = 0; j < nc; ++j)
      mtx[i][j] = 0;
    mtx[i][i] = 1;
  }
  return mtx;
}
void freeMatrix(double **mtx){
  if(mtx == NULL) return; //nothing to free...
  void *indi = (void*)mtx - sizeof(int);
  int nr = ((int*)indi)[0];
  if(nr) for (int i = 0; i < nr; ++i) free(mtx[i]);
  else free(mtx[0]);
  free(indi);
}

//simple operations
void addMatrix(double **m1, double** m2, int n, int m, double **out){
  for (int i = 0; i < n; ++i)
    sumVectors(m1[i], m2[i], out[i], m);
}
void substractMatrix(double **m1, double** m2, int n, int m, double **out){
  for (int i = 0; i < n; ++i)
    substractVectors(m1[i], m2[i], out[i], m);
}
//multiply nxm with mxp
void multiplyMatrix(double **m1, double **m2, int n, int m, int p, double **out){
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < m; ++j)
    {
      double c = 0;
      for (int k = 0; k < p; ++k)
      {
        c += m1[i][k] * m2[k][j];
      }
      out[i][j] = c;
    }
  }
}
void multiplyMatrixTransposed(double **m1, double **m2, int n, int m, int p, double **out){
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < m; ++j)
    {
      double c = 0;
      for (int k = 0; k < p; ++k)
      {
        c += m1[i][k] * m2[j][k];
      }
      out[i][j] = c;
    }
  }
}
//vector matrix
void multiplyMatrixVector(double **mat, double *v, int n, int m, double *out){
  for (int i = 0; i < n; i++) {
    out[i] = dotproduct(mat[i], v, m);
  }
}