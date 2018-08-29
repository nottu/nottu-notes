#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "optimization.h"

double to_degrees(double radians) {
  return radians * (180.0/M_PI);
}

int main(int argc, char** argv) {
  double xi[] = {746,  629, 1571, 155};
  double yi[] = {1393, 375, 259,  987};
  double ti[] = {-4, 15, 52};
  double si[] = {8, 4.6, 3.3, 4.0};
  double d4   = 864.3;

  double *xy = newVector(2);
  xy[0] = 780; xy[1] = 600;

  double func(double *v, int n){
    double f = 0;

    double degerees = 90 + to_degrees(atan2( (v[1] - yi[0]), (xi[0] - v[0]) ));
    double tmp = ti[0] - degerees;
    f += SQUARE(tmp/si[0]);

    degerees = 90 - to_degrees(atan2((v[1] - yi[1]), -(xi[1] - v[0]) ));
    tmp = ti[1] - degerees;
    f += SQUARE(tmp/si[1]);

    degerees = 90 + to_degrees(atan2((v[1] - yi[2]), (xi[2] - v[0]) ));
    tmp = ti[2] - degerees;
    f += SQUARE(tmp/si[2]);

    tmp = SQUARE(v[1] - yi[3]) + SQUARE(v[0] - xi[3]) - SQUARE(d4);
    f += SQUARE(tmp/si[3]);

    return f;
  }
  double *f_grad(double *x, int n, double *grad){
    return finete_difference(x, func, n, 0.0001, grad);
  }

  FuncInfo info;
  info.function = func;
  info.gradient = f_grad;

  double toler = 1E-5;

  printf("Ini\t%g\n", func(xy, 2));
////
  double err = optimize_function(info, StepBacktrack, xy, 2, 100000, toler, toler, toler, 1);

  printf("Error\t%g\n", err);
  printf("x(");printVector(xy, 2);printf(")\n");
  printf("\n");

  printf("%g\n",  90 + to_degrees(atan2( (xy[1] - yi[0]), (xi[0] - xy[0]) )));
  printf("%g\n",  90 + to_degrees(atan2((xy[1] - yi[1]), -(xi[1] - xy[0]) )) );
  printf("%g\n",  90 - to_degrees(atan2((xy[1] - yi[2]), (xi[2] - xy[0]) )) );
  printf("%g\n",  pow(SQUARE(xy[1] - yi[2])+ SQUARE(xi[2] - xy[0]),1.0/2.0));
//
  free(xy);

  return 0;
}