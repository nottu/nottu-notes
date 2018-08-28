#include <stdio.h>
#include <stdlib.h>

#include "optimization.h"
#include "ej2_funcs.h"


void optimizeGaussianAdjust(int ***hist, double *alphas, double **mus, int x, int y, int z, double sigma, int n){
  double fx = gaussianAdjust(hist, x, y, z, alphas, mus, sigma, n);
  printf("FX %g\n", fx);
  double *g = gaussianAdjust_alpha_gradient(hist, x, y, z, alphas, mus, sigma, n, newVector(n));
  printf("Gradient\t");printVector(g, n);
  double **h = gaussianAdjust_alpha_hessian(hist, x, y, z, alphas, mus, sigma, n, allocMatrix(n, n));

  int max_iter = 100;

  for (int i = 0; i < 1000; ++i)
  {
    printf("\n\n\nalp\n");
    doglegOptimize_alpha(hist, x, y, z, alphas, mus, sigma, n, max_iter, 1E-4, 10);

    // normalizeVector(alphas, n);
    printf("\n\n\nCH1\n");
    doglegOptimize_mu(hist, x, y, z, alphas, mus, sigma, n, max_iter, 1E-4, 10, 0);

    printf("\n\n\nCH2\n");
    doglegOptimize_mu(hist, x, y, z, alphas, mus, sigma, n, max_iter, 1E-4, 10, 1);

    printf("\n\n\nCH3\n");
    doglegOptimize_mu(hist, x, y, z, alphas, mus, sigma, n, max_iter, 1E-6, 10, 2);
  }
  gaussianAdjust_alpha_gradient(hist, x, y, z, alphas, mus, sigma, n, g);
  printf("Gradient\t");printVector(g, n);
  fx = gaussianAdjust(hist, x, y, z, alphas, mus, sigma, n);
  printf("FX %g\n", fx);
  printf("alphas\t");printVector(alphas, n);
  printf("MUS\n");
  printMatrix(mus, n, 3);

  freeMatrix(h);
  free(g);
}

void testGaussianAdjust(double *alphas1, double *alphas2, double **mus1, double **mus2, int x, int y, int z, double sigma, int n){
  //hardcode image file
  FILE *f = fopen("segmentation.txt", "w");
  fprintf(f, "%i %i %i\n", x, y, z);
  double *cvec = newVector(DIM); //[r, g, b]
  double *c_mu = newVector(DIM);
  double eps = 0.01;
  for (int i = 0; i < x; ++i){
  cvec[0] = i;// * 255/(x-1);
    for (int j = 0; j < y; ++j){
      cvec[1] = j;// * 255/(x-1);
      for (int k = 0; k < z; ++k){
        cvec[2] = j;
        // evalGaussianMix(alphas, mus, cvec, sigma, n);
        double f1 = evalGaussianMix(alphas1, mus1, cvec, sigma, n);
        double f2 = evalGaussianMix(alphas2, mus2, cvec, sigma, n);
        double F1 = (f1 + eps) / (f1 + f2 + 2 * eps);
        double F2 = (f2 + eps) / (f1 + f2 + 2 * eps);
        fprintf(f, "%i\n", F1 > F2 ? 1 : 2);
      }
    }
  }
  fclose(f);
  free(cvec);
  free(c_mu);
}

void doGaussianAdjust(int ***hist1, int ***hist2, double *alphas1, double *alphas2, double **mus1, double **mus2, int x, int y, int z, double sigma, int n){
  optimizeGaussianAdjust(hist1, alphas1, mus1, x, y, z, sigma, n);
  optimizeGaussianAdjust(hist2, alphas2, mus2, x, y, z, sigma, n);

  FILE *f = fopen("out.txt", "w");
  if(f == NULL) return;
  fprintf(f, "ALPHAS_1\n");
  for (int i = 0; i < n; ++i) fprintf(f, "%g ", alphas1[i]);
  fprintf(f, "\nMUS_1\n");
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < 3; ++j) {
      fprintf(f, "%g ", mus1[i][j]); 
    }
    fprintf(f, "\n");
  }
  fprintf(f, "ALPHAS_2\n");
  for (int i = 0; i < n; ++i) fprintf(f, "%g ", alphas1[i]);
  fprintf(f, "\nMUS_2\n");
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < 3; ++j) {
      fprintf(f, "%g ", mus1[i][j]); 
    }
    fprintf(f, "\n");
  }
  testGaussianAdjust(alphas1, alphas2, mus1, mus2, x, y, z, sigma, n);
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
  for (int i = 0; i < n; ++i) resetVector(mus1[i], 3, (i+1));

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