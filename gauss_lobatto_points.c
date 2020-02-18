#include <math.h>

void gauss_lobatto_points(double *xi, int nt){
   int i;
   for(i=0; i<=nt; i++){
      xi[i] = -cos(M_PI/nt*i);
   }
}
