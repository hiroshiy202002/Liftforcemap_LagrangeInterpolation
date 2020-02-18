#include <stdio.h>
#include <math.h>
#include "lagrange_interpolation_tools.h"

// #define YMAX (1.0-(1.0+8.0/20)/6.7)
#define YMAX 0.791045
#define N 160
#define NT 14
#define ND 2

void inputdata(double*, double*, double*, int);
void lagint(double*, double*, double*, double*, double*, double*, int, int, int);
void outputdata(double*, double*, double*, int, double);

int main(void){
   double Fx[(NT+1)*(NT+1)], Fy[(NT+1)*(NT+1)], Fz[(NT+1)*(NT+1)];
   double IntpltdFx[(N+1)*(N+1)];
   double IntpltdFr[(N+1)*(N+1)];
   double IntpltdFt[(N+1)*(N+1)];
   inputdata(Fx, Fy, Fz, NT);
   lagint(Fx, Fy, Fz, IntpltdFx, IntpltdFr, IntpltdFt, ND, NT, N);
   outputdata(IntpltdFx, IntpltdFr, IntpltdFt, N, YMAX);
   return(0);
}

void inputdata(double *fx, double *fy, double *fz, int nt){
   int ntp1 = nt+1;
   int jt, kt;
   FILE *fp;
   char inputfilename[] = "lat_at_nodes.dat";
   fp = fopen(inputfilename, "r");
      for(jt=nt/2; jt<=nt; jt++){
         for(kt=nt/2; kt<=jt; kt++){
            fscanf(fp, "%lf %lf %lf",
               &fx[jt+kt*ntp1], &fy[jt+kt*ntp1], &fz[jt+kt*ntp1]);
         }
      }
   fclose(fp);
   for(jt=nt/2+1; jt<=nt; jt++){
      for(kt=nt/2; kt<=jt; kt++){
         fx[kt+jt*ntp1] = fx[jt+kt*ntp1];
         fz[kt+jt*ntp1] = fy[jt+kt*ntp1];
         fy[kt+jt*ntp1] = fz[jt+kt*ntp1];
      }
   }
   for(kt=nt/2; kt<=nt; kt++){
      for(jt=nt/2; jt<=nt; jt++){
         fx[(nt-jt)+kt*ntp1] = fx[jt+kt*ntp1];
         fy[(nt-jt)+kt*ntp1] = -fy[jt+kt*ntp1];
         fz[(nt-jt)+kt*ntp1] = fz[jt+kt*ntp1];
         fx[jt+(nt-kt)*ntp1] = fx[jt+kt*ntp1];
         fy[jt+(nt-kt)*ntp1] = fy[jt+kt*ntp1];
         fz[jt+(nt-kt)*ntp1] = -fz[jt+kt*ntp1];
         fx[(nt-jt)+(nt-kt)*ntp1] = fx[jt+kt*ntp1];
         fy[(nt-jt)+(nt-kt)*ntp1] = -fy[jt+kt*ntp1];
         fz[(nt-jt)+(nt-kt)*ntp1] = -fz[jt+kt*ntp1];
      }
   }
   for(kt=0; kt<=nt; kt++) fy[nt/2+kt*ntp1] = 0.0;
   for(jt=0; jt<=nt; jt++) fz[jt+nt/2*ntp1] = 0.0;
}

void lagint(double *fx, double *fy, double *fz,
            double *intpltdfx, double *intpltdfr, double *intpltdft,
            int nd, int nt, int n){
   int ntp1 = nt+1, np1 = n+1;
   double xi[ntp1];
   double intpltdfy[np1*np1], intpltdfz[np1*np1];
   double y, z, hjthkt, r;
   int j, k, jt, kt;
   int i;
   gauss_lobatto_points(xi, nt);
   for(k=0; k<=n; k++){
      z = -1.0+2.0/N*k;
      for(j=0; j<=n; j++){
         y = -1.0+2.0/N*j;
         intpltdfx[j+k*np1] = 0.0;
         intpltdfy[j+k*np1] = 0.0;
         intpltdfz[j+k*np1] = 0.0;
         for(kt=0; kt<=nt; kt++){
            for(jt=0; jt<=nt; jt++){
               hjthkt = hk(0, jt, y, xi, nd, nt)*hk(0, kt, z, xi, nd, nt);
               intpltdfx[j+k*np1] += fx[jt+kt*ntp1]*hjthkt;
               intpltdfy[j+k*np1] += fy[jt+kt*ntp1]*hjthkt;
               intpltdfz[j+k*np1] += fz[jt+kt*ntp1]*hjthkt;
            }
         }
         if(j==n/2 && k==n/2){
            intpltdfr[j+k*np1] = 0.0;
            intpltdft[j+k*np1] = 0.0;
            continue;
         }
         r = sqrt(y*y+z*z);
         intpltdfr[j+k*np1] = y/r*intpltdfy[j+k*np1]+z/r*intpltdfz[j+k*np1];
         intpltdft[j+k*np1] = -z/r*intpltdfy[j+k*np1]+y/r*intpltdfz[j+k*np1];
      }
      printf("%d/%d\n", k, n);
   }
}

void outputdata(double *intpltdfx, double *intpltdfr, double *intpltdft,
                                                         int n, double ymax){
   int np1 = n+1;
   double y, z;
   int j, k;
   FILE *fpfx, *fpfr, *fpft;
   char inputfilename1[] = "IntpltdFx.dat";
   char inputfilename2[] = "IntpltdFr.dat";
   char inputfilename3[] = "IntpltdFt.dat";
   fpfx = fopen(inputfilename1, "w");
   fpfr = fopen(inputfilename2, "w");
   fpft = fopen(inputfilename3, "w");
      fprintf(fpfx, "# y z Fx\n");
      fprintf(fpfr, "# y z Fr\n");
      fprintf(fpft, "# y z Ft\n");
      for(k=0; k<=n; k++){
         z = -1.0+2.0/N*k;
         for(j=0; j<=n; j++){
            y = -1.0+2.0/N*j;
            fprintf(fpfx, "%19.16e %19.16e %19.16e\n",
                              y*ymax, z*ymax, intpltdfx[j+k*np1]);
            fprintf(fpfr, "%19.16e %19.16e %19.16e\n",
                              y*ymax, z*ymax, intpltdfr[j+k*np1]);
            fprintf(fpft, "%19.16e %19.16e %19.16e\n",
                              y*ymax, z*ymax, intpltdft[j+k*np1]);
         }
         fprintf(fpfx, "\n");
         fprintf(fpfr, "\n");
         fprintf(fpft, "\n");
      }
   fclose(fpfx);
   fclose(fpfr);
   fclose(fpft);
}
