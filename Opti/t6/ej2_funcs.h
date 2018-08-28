#ifndef EJ2_FUNCS
#define EJ2_FUNCS
#define DIM 3
/* Histogram reading stuff */
int ***read3DHist(char* filename, int *x, int *y, int *z);
void free3DHist(int ***hist_3d, int xv, int yv, int zv);
void print3Dhist(int ***hist, int x, int y, int z);

/* Eval Functions*/
double evalGaussianMix(double *alphas, double **mus, double *c, double sigma, int n);
double gaussianAdjust(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n);

/* Gradient Functions*/
double* gaussianAdjust_alpha_gradient(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n, double *g);
double* gaussianAdjust_mu_gradient(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n, double *g, int channel);

/* Hessian Functions*/
double **gaussianAdjust_alpha_hessian(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n, double **h);
double **gaussianAdjust_mu_hessian(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n, double **h, int channel);

/*Dog Leg Alpha stuff*/
double doglegOptimize_alpha(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n, int max_iter, double tg, double reg_szM);
double doglegOptimize_mu(int ***hist, int x, int y, int z, double *alphas, double **mus, double sigma, int n, int max_iter, double tg, double reg_szM, int channel);
#endif