#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ej2_funcs.h"
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
  double *c_mu = newVector(DIM);
  for (int i = 0; i < n; ++i){
    substractVectors(c, mus[i], c_mu, DIM);
    double exp_term = -1 * SQUARE(vectorNorm(c_mu, DIM, 2));
    exp_term /= 2 * SQUARE(sigma);
    val+= alphas[i] * exp(exp_term);
  }
  free(c_mu);
  return val;
}

double gaussianAdjust(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n){
  double val = 0;
  double *cvec = newVector(DIM); //[r, g, b]
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
  free(cvec);
  return val;
}

//* ALPHA STUFF*//
double* gaussianAdjust_alpha_gradient(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n, double *gradient){
  // double *gradient = newVector(n);
  double *cvec = newCleanVector(3); //[r, g, b]
  double *c_mu = newCleanVector(3);
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
          substractVectors(cvec, mus[l], c_mu, 3);
          double e1  = vectorNorm(c_mu, 3, 2);
                 e1 *= e1 * -1;
                 e1 /=  2 * SQUARE(sigma);
          val += -2 * (hc - gm) * exp(e1);
        }
      }
    }
    gradient[l] = val;
  }
  free(c_mu);
  free(cvec);
  return gradient;
}

double **gaussianAdjust_alpha_hessian(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n, double **hessian){
  // double **hessian = allocMatrix(n, n);
  double *cvec = newCleanVector(3); //[r, g, b]
  double *c_mu = newCleanVector(3);
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
            substractVectors(cvec, mus[l], c_mu, 3);
            double e1 =  exp(-1 * SQUARE(vectorNorm(c_mu, DIM, 2))  / (2*SQUARE(sigma)) );
            substractVectors(cvec, mus[m], c_mu, 3);
            val += e1 * exp(-1 * SQUARE(vectorNorm(c_mu, DIM, 2)) / (2*SQUARE(sigma)) );
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

double doglegOptimize_alpha(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n, int max_iter, double tg, double reg_szM){
  double *g = gaussianAdjust_alpha_gradient(hist, x, y, z, alphas, mus, sigma, n, newVector(n));
  double **H = gaussianAdjust_alpha_hessian(hist, x, y, z, alphas, mus, sigma, n, allocMatrix(n, n));

  double fx;// = gaussianAdjust(hist, x, y, z, alphas, mus, sigma, n);
  double norm2 = vectorNorm(g, n, 2);
  double reg_sz = reg_szM;
  int iter = 0;

  double *pu = newVector(n);
  double *dir = newVector(n);
  double *alphas1 = newVector(n);
  double phi, taylorev;
  while(norm2 >= tg && iter < max_iter){
    fx = gaussianAdjust(hist, x, y, z, alphas, mus, sigma, n);
    double a = cauchy_point(n, g, H);
    copyVector(g, pu, n);
    scaleVector(pu, n, a);
    double pu_norm = vectorNorm(pu, n ,2);
    copyVector(pu, dir, n);
    if(pu_norm < reg_sz){
      if(isdefpos(H, n)) {
        //do dogleg direction
        // printf("DOGLEG!\n");
        doglegDirection(n, g, H, pu, reg_sz, dir);
      }
    } else {
      //make smaller
      scaleVector(dir, n, reg_sz/pu_norm);
    }
    sumVectors(alphas, dir, alphas1, n); //x = x + dir
    taylorev = fx - taylorEval(fx, n, g, H, dir);
    double fx1 = gaussianAdjust(hist, x, y, z, alphas1, mus, sigma, n);
    phi = (fx - fx1 ) / (taylorev);
    //....
    // printf("Iter %i:%i \tf(x): %g\t||g||: %g\t reg_sz %lg\t phi: %g\t taylor %g\n", iter, max_iter, fx, norm2, reg_sz, phi, taylorev);
    iter++;
    if(isnan(phi)) break;
    if(phi < 0.25) {
      reg_sz /= 4; //bad model
      if(reg_sz < 1E-6) break;
      // printf("\tBAD MODEL\tphi %g\n", phi);
    }
    if(phi < 0) continue;
    double dirnorm = vectorNorm(dir, n, 2);
    if(phi >= 0.75 && dirnorm + 0.001 > reg_sz ){
      reg_sz = 2*reg_sz < reg_szM ? 2*reg_sz : reg_szM;
    }
    copyVector(alphas1, alphas, n);
    g = gaussianAdjust_alpha_gradient(hist, x, y, z, alphas, mus, sigma, n, g);
    H = gaussianAdjust_alpha_hessian(hist, x, y, z, alphas, mus, sigma, n, H);
    norm2 = vectorNorm(g, n, 2);
  }
  printf("Exit Iter %i:%i \tf(x): %g\t||g||: %g\t reg_sz %lg\t phi: %g\t taylor %g\n", iter, max_iter, fx, norm2, reg_sz, phi, taylorev);

  free(g);
  free(alphas1);
  free(pu);
  free(dir);
  freeMatrix(H);
  return fx;
}

//* MU STUFF*//
double* gaussianAdjust_mu_gradient(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n, double *gradient, int channel){
  double *cvec = newVector(DIM); //[r, g, b]
  double *c_mu = newVector(DIM);
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

          substractVectors(cvec, mus[l], c_mu, 3);
          val += -1 * (hc - gm) * alphas[l] * exp(-1 * SQUARE(vectorNorm(c_mu, DIM, 2)) / (2*SQUARE(sigma)) ) * (c_mu[channel] - mus[l][channel]) / SQUARE(sigma);
        }
      }
    }
    gradient[l] = val;
  }
  free(c_mu);
  free(cvec);
  return gradient;
}


double **gaussianAdjust_mu_hessian(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n, double **hessian, int channel){
  // double **hessian = allocMatrix(n, n);
  double *cvec = newVector(DIM); //[r, g, b]
  double *c_mu = newVector(DIM);
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
            substractVectors(cvec, mus[l], c_mu, 3);
            double e1 =  exp(-1 * SQUARE(vectorNorm(c_mu, DIM, 2)) / (2*SQUARE(sigma)) ) * (c_mu[channel] - mus[l][channel]) / SQUARE(sigma);
            substractVectors(cvec, mus[m], c_mu, 3);
            val += e1 * exp(-1 * SQUARE(vectorNorm(c_mu, DIM, 2)) / (2*SQUARE(sigma)) ) * (c_mu[channel] - mus[m][channel]) / SQUARE(sigma);
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

double doglegOptimize_mu(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n, int max_iter, double tg, double reg_szM, int channel){
  double *g = gaussianAdjust_mu_gradient(hist, x, y, z, alphas, mus, sigma, n, newVector(n), channel);
  double **H = gaussianAdjust_mu_hessian(hist, x, y, z, alphas, mus, sigma, n, allocMatrix(n, n), channel);

  double fx;// = gaussianAdjust(hist, x, y, z, alphas, mus, sigma, n);
  double norm2 = vectorNorm(g, n, 2);
  double reg_sz = reg_szM;
  int iter = 0;

  double *pu = newVector(n);
  double *dir = newVector(n);
  double **mus1 = allocMatrix(n, 3);
  double phi, taylorev;
  while(norm2 >= tg && iter < max_iter){
    copyMatrix(mus, mus1, n, 3);
    fx = gaussianAdjust(hist, x, y, z, alphas, mus, sigma, n);
    iter++;
    double a = cauchy_point(n, g, H);
    copyVector(g, pu, n);
    scaleVector(pu, n, a);
    double pu_norm = vectorNorm(pu, n ,2);
    copyVector(pu, dir, n);
    if(pu_norm < reg_sz){
      if(isdefpos(H, n)) {
        //do dogleg direction
        // printf("DOGLEG!\n");
        doglegDirection(n, g, H, pu, reg_sz, dir);
      }
    } else {
      //make smaller
      scaleVector(dir, n, reg_sz/pu_norm);
    }
    // sumVectors(alphas, dir, alphas1, n); //x = x + dir
    for (int i = 0; i < n; ++i){
      mus1[i][channel] = mus[i][channel] + dir[i];
    }
    taylorev = fx - taylorEval(fx, n, g, H, dir);
    double fx1 = gaussianAdjust(hist, x, y, z, alphas, mus1, sigma, n);
    phi = (fx - fx1 ) / (taylorev);
    //....
    // printf("Iter %i:%i \tf(x): %g\t||g||: %g\t reg_sz %lg\t phi: %g\t taylor %g\n", iter, max_iter, fx, norm2, reg_sz, phi, taylorev);
    if(phi < 0.25) {
      reg_sz /= 4; //bad model
      if(reg_sz < 1E-6) break;
      // printf("\tBAD MODEL\tphi %g\n", phi);
    }
    if(phi < 0) continue;
    double dirnorm = vectorNorm(dir, n, 2);
    if(phi >= 0.75 && dirnorm + 0.001 > reg_sz ){
      reg_sz = 2*reg_sz < reg_szM ? 2*reg_sz : reg_szM;
    }
    copyMatrix(mus1, mus, n, 3);
    g = gaussianAdjust_mu_gradient(hist, x, y, z, alphas, mus, sigma, n, g, channel);
    H = gaussianAdjust_mu_hessian(hist, x, y, z, alphas, mus, sigma, n, H, channel);
    norm2 = vectorNorm(g, n, 2);
  }
  printf("Exit Iter %i:%i \tf(x): %g\t||g||: %g\t reg_sz %lg\t phi: %g\t taylor %g\n", iter, max_iter, fx, norm2, reg_sz, phi, taylorev);

  free(g);
  freeMatrix(mus1);
  free(pu);
  free(dir);
  freeMatrix(H);
  return fx;
}
