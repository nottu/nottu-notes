#ifndef OPTI_UTILS
#define OPTI_UTILS
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

double *readVector(char* name, int* sz);
double **readMtx(char* name, int* nr, int* nc);
#endif