void chebyshev_polynomials(double *Tn, double x, int nd, int nt){
   int ndp1 = nd+1;
   int id, it;
   for(id=0; id<=nd; id++){
      Tn[id+0*ndp1] = 0.0;
      Tn[id+1*ndp1] = 0.0;
   }
   Tn[0+0*ndp1] = 1.0;
   Tn[0+1*ndp1] = x;
   Tn[1+1*ndp1] = 1.0;
   for(it=2; it<=nt; it++){
      Tn[0+it*ndp1] = 2.0*x*Tn[0+(it-1)*ndp1] - Tn[0+(it-2)*ndp1];
      for(id=1; id<=nd; id++){
         Tn[id+it*ndp1] = 2.0*id*Tn[(id-1)+(it-1)*ndp1]
                        + 2.0*x*Tn[id+(it-1)*ndp1] - Tn[id+(it-2)*ndp1];
      }
   }
}
