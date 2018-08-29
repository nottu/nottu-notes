#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv){
  FILE *f = fopen("convex.txt", "w");
  int n = 1000;
  fprintf(f, "%i\n", n);
  for (int i = 1; i < n + 1; ++i)
  {
    fprintf(f, "%0.10lf\n", 1.0);
  }
  fclose(f);
}