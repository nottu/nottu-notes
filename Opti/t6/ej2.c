#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "optimization.h"

int ***read3DHist(char* filename, int *x, int *y, int *z){
  FILE *f = fopen(filename, "r");
  if(f != NULL){
    int val = fscanf(f, "%i %i %i", x, y, z);
    int xv = *x, yv = *y, zv = *z;
    int ***hist_3d = (int***)malloc(sizeof(int **) * xv);
    for (int i = 0; i < xv; ++i)
    {
      hist_3d[i] = (int**)malloc(sizeof(int*) * yv);
      for (int j = 0; j < yv; ++j)
      {
        hist_3d[i][j] = (int*)malloc(sizeof(int) * zv);
        for (int k = 0; k < zv; ++k)
        {
          fscanf(f, "%d", &hist_3d[i][j][k]);
        }
      }
    }
    fclose(f);
    return hist_3d;
  }
  exit(1);
}
void free3DHist(int ***hist_3d, int xv, int yv, int zv){
  for (int i = 0; i < xv; ++i) {
    for (int j = 0; j < yv; ++j) {
      free(hist_3d[i][j]);
    }
    free(hist_3d[i]);
  }
  free(hist_3d);
}
void print3Dhist(int ***hist, int x, int y, int z){
  for (int i = 0; i < x; ++i){
    for (int j = 0; j < y; ++j){
      for (int k = 0; k < z; ++k){
        printf("%d\t", hist[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }
}

double evalGaussianMix(double *alphas, double **mus, double *c, double sigma, int n){
  double val = 0;
  double *c_mu = newVector(3);
  for (int i = 0; i < n; ++i){
    substractVectors(c, mus[i], c_mu, 3);
    double exp_term = -1 * SQUARE(vectorNorm(c_mu, 3, 2));
    exp_term /= 2 * SQUARE(sigma);
    val+= alphas[i] * exp(exp_term);
  }
  free(c_mu);
  return val;
}

double gaussianAdjust(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n){
  double val = 0;
  double *cvec = newVector(3); //[r, g, b]
  for (int i = 0; i < x; ++i){
    cvec[0] = i;
    for (int j = 0; j < y; ++j){
      cvec[1] = j;
      for (int k = 0; k < z; ++k){
        cvec[2] = k;
        double gm = evalGaussianMix(alphas, mus, cvec, sigma, n);
        double hc = hist[i][j][k];
        val += SQUARE(hc - gm);
      }
    }
  }
  return val;
}

double* gaussianAdjust_alpha_gradient(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n){
  double *gradient = newVector(n);
  double *cvec = newVector(3); //[r, g, b]
  double *c_mu = newVector(3);
  for (int l = 0; l < n; ++l){

    double val = 0;
    // create c
    for (int i = 0; i < x; ++i){
      cvec[0] = i;
      for (int j = 0; j < y; ++j){
        cvec[1] = j;
        for (int k = 0; k < z; ++k){
          cvec[2] = k;
          double gm = evalGaussianMix(alphas, mus, cvec, sigma, n);
          double hc = hist[i][j][k];

          substractVectors(cvec, mus[l], c_mu, n);
          val += -2 * (hc - gm) * exp(-1 * SQUARE(vectorNorm(c_mu, 3, 2)));
        }
      }
    }
    gradient[l] = val;
  }
  free(c_mu);
  free(cvec);
  return gradient;
}

double **gaussianAdjust_alpha_hessian(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n){
  double **hessian = allocMatrix(n, n);
  double *cvec = newVector(3); //[r, g, b]
  double *c_mu = newVector(3);
  for (int l = 0; l < n; ++l){
    for (int m = l; m < n; ++m){

      double val = 0;
      // create c
      for (int i = 0; i < x; ++i){
        cvec[0] = i;
        for (int j = 0; j < y; ++j){
          cvec[1] = j;
          for (int k = 0; k < z; ++k){
            cvec[2] = k;
            substractVectors(cvec, mus[l], c_mu, n);
            double e1 =  exp(-1 * SQUARE(vectorNorm(c_mu, 3, 2)));
            substractVectors(cvec, mus[m], c_mu, n);
            val += e1 * exp(-1 * SQUARE(vectorNorm(c_mu, 3, 2)));
          }
        }
      }

      hessian[l][m] = val;
      hessian[m][l] = val; //hessian is symmetrical
    }
  }
  free(c_mu);
  free(cvec);
  return hessian;
}


void optimizeGaussianAdjust(int ***hist, double *alphas, double **mus, int x, int y, int z, double sigma, int n){
  double fx = gaussianAdjust(hist, x, y, z, alphas, mus, sigma, n);
  printf("%g\n", fx);
  double *g = gaussianAdjust_alpha_gradient(hist, x, y, z, alphas, mus, sigma, n);
  printf("gradient:\t"); printVector(g, n);

  double **h = gaussianAdjust_alpha_hessian(hist, x, y, z, alphas, mus, sigma, n);
  printf("hessian\n"); printMatrix(h, n, n);

  freeMatrix(h);
  free(g);
}

void doGaussianAdjust(int ***hist1, int ***hist2, double *alphas1, double *alphas2, double **mus1, double **mus2, int x, int y, int z, double sigma, int n){
  optimizeGaussianAdjust(hist1, alphas1, mus1, x, y, z, sigma, n);
  optimizeGaussianAdjust(hist2, alphas2, mus2, x, y, z, sigma, n);
}

int main(int argc, char **argv){
  if(argc < 5) return(1);
  int n = atoi(argv[1]);
  int sigma = atoi(argv[2]);
  int x, y, z;
  int ***hist1 =  read3DHist(argv[3], &x, &y, &z);
  int ***hist2 =  read3DHist(argv[4], &x, &y, &z);
  // print3Dhist(hist1, x, y, z);

  double *alphas1 = filledVector(n, 1);
  normalizeVector(alphas1, n);
  double **mus1 = allocMatrix(n, 3);
  for (int i = 0; i < n; ++i) resetVector(mus1[i], 3, i);

  double *alphas2 = filledVector(n, 1);
  normalizeVector(alphas2, n);
  double **mus2 = allocMatrix(n, 3);
  for (int i = 0; i < n; ++i) resetVector(mus2[i], 3, i);

  doGaussianAdjust(hist1, hist2, alphas1, alphas2, mus1, mus2, x, y, z, sigma, n);

  free(alphas1);
  free(alphas2);

  freeMatrix(mus1);
  freeMatrix(mus2);

  free3DHist(hist1, x, y, z);
  free3DHist(hist2, x, y, z);
  return 0;
}