#include "opti_utils.h"

double *readVector(char* name, int* sz){
  FILE *f = fopen(name, "rb");
  if (!f) return NULL;
  fread(sz, sizeof(int), 1, f);
  double *vect = newVector(*sz);
  for (int i = 0; i < *sz; ++i) {
    fread(vect, sizeof(double), *sz, f);
  }
  fclose(f);
  return vect;
}
double **readMtx(char* name, int* nr, int* nc){
  FILE *f = fopen(name, "rb");
  if (!f) return NULL;
  fread(nr, sizeof(int), 1, f);
  fread(nc, sizeof(int), 1, f);

  double **mtx = allocMatrix(*nr, *nc);
  for (int i = 0; i < *nr; ++i) {
    fread(mtx[i], sizeof(double), (unsigned int)*nc, f);
  }
  fclose(f);
  return mtx;
}