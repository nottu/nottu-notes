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

  for (int i = 0; i < max_iter; ++i)
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

void doGaussianAdjust(int ***hist1, int ***hist2, double *alphas1, double *alphas2, double **mus1, double **mus2, int x, int y, int z, double sigma, int n){
  optimizeGaussianAdjust(hist1, alphas1, mus1, x, y, z, sigma, n);
  // optimizeGaussianAdjust(hist2, alphas2, mus2, x, y, z, sigma, n);
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
  // alphas1[0] = 133444;
  // alphas1[1] = -391156;
  normalizeVector(alphas1, n);
  double **mus1 = allocMatrix(n, 3);
  for (int i = 0; i < n; ++i) resetVector(mus1[i], 3, (i+1));
  // for (int i = 0; i < n; ++i) {
    // mus1[0][0] =  1;
    // mus1[0][1] = -1;
    // mus1[0][2] =  3.3;

    // mus1[1][0] =  1;
    // mus1[1][1] = -2;
    // mus1[1][2] =  3.6;
    
  // }

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