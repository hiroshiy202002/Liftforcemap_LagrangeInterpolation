void chebyshev_polynomials(double*, double, int, int);

double hk(int id, int it, double x, double *xi, int nd, int nt){
   double hk;
   double tn1[(nd+1)*(nt+1)], tn2[(nd+1)*(nt+1)];
   double c[nt+1];
   int ndp1 = nd+1;
   double tmp;
   int i;
   for(i=1; i<nt; i++){
      c[i] = 1.0;
   }
   c[0] = 2.0;
   c[nt] = 2.0;
   chebyshev_polynomials(tn1, xi[it], nd, nt);
   chebyshev_polynomials(tn2, x, nd, nt);
   tmp = 0.0;
   for(i=0; i<=nt; i++){
      tmp += tn1[0+i*ndp1]*tn2[id+i*ndp1]/(c[it]*c[i]);
   }
   hk = 2.0/nt*tmp;
   return hk;
}
