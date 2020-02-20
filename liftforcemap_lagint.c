#include <stdio.h>
#include <math.h>
#include "lagrange_interpolation_tools.h"

// #define YMAX (1.0-(1.0+8.0/20)/6.7)
#define YMAX 0.791045
#define N 160
#define NT 14
#define ND 2

#define LARGEVAL 1.0e7
#define MYFMT "%19.16e"

void inputdata(double*, double*, double*, int);
void lagint(double*, double*, double*, double*, double*, double*, int, int, int);
void search_equilibriumpoints(double*, double*, int, int, double*, double*);
void equilibriumpoint_by_newton_raphson_method
            (double*, double*, int, int, double*, double*);
void outputdata(double*, double*, double*, int, double*, double*, double);

int main(void){
   double Fx[(NT+1)*(NT+1)], Fy[(NT+1)*(NT+1)], Fz[(NT+1)*(NT+1)];
   double IntpltdFx[(N+1)*(N+1)];
   double IntpltdFr[(N+1)*(N+1)];
   double IntpltdFt[(N+1)*(N+1)];
   double Ye[8], Ze[8];
   int i;
   inputdata(Fx, Fy, Fz, NT);
   lagint(Fx, Fy, Fz, IntpltdFx, IntpltdFr, IntpltdFt, ND, NT, N);
   search_equilibriumpoints(Fy, Fz, ND, NT, Ye, Ze);
   outputdata(IntpltdFx, IntpltdFr, IntpltdFt, N, Ye, Ze, YMAX);
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

void search_equilibriumpoints(double *fy, double *fz, int nd, int nt,
                                                double *ye, double *ze){
   double y, z;
   int counter=0, i, j, k, flag;
   int n = 15;
   double eps = 1.0e-7;
   for(k=0; k<=n; k++){
      z = 1.0/n*k;
      for(j=1; j<=n; j++){
         y = 1.0/n*j;
         equilibriumpoint_by_newton_raphson_method(fy, fz, nd, nt, &y, &z);
         if(y==LARGEVAL) continue;
         if(y<eps && z<eps) continue;
         // printf("y="MYFMT" z="MYFMT"\n", y, z);
         if(counter==0){
            ye[counter] = y;
            ze[counter] = z;
            printf("ye%d=%lf ze%d=%lf\n",
               counter, ye[counter], counter, ze[counter]);
            counter += 1;
         }else{
            flag = 0;
            for(i=0; i<counter; i++){
               if(fabs(ye[i]-y)<eps && fabs(ze[i]-z)<eps) flag = 1;
            }
            if(flag==0){
               ye[counter] = y;
               ze[counter] = z;
               printf("ye%d=%lf ze%d=%lf\n",
                  counter, ye[counter], counter, ze[counter]);
               counter += 1;
            }
         }
      }
   }
   for(i=counter; i<8; i++){
      ye[i] = LARGEVAL;
      ze[i] = LARGEVAL;
   }
}

void equilibriumpoint_by_newton_raphson_method
            (double *fy, double *fz, int nd, int nt, double *y, double *z){
   int ntp1 = nt+1;
   double xi[ntp1];
   double hjthkt, ddy_hjthkt, ddz_hjthkt;
   double Fy, Fz, ddy_Fy, ddy_Fz, ddz_Fy, ddz_Fz, J;
   int counter, jt, kt;
   double eps = 1.0e-30, tmp;
   gauss_lobatto_points(xi, nt);
   for(counter=1;;counter++){
      Fy = 0.0;
      Fz = 0.0;
      ddy_Fy = 0.0;
      ddy_Fz = 0.0;
      ddz_Fy = 0.0;
      ddz_Fz = 0.0;
      for(kt=0; kt<=nt; kt++){
         for(jt=0; jt<=nt; jt++){
            hjthkt = hk(0, jt, *y, xi, nd, nt)*hk(0, kt, *z, xi, nd, nt);
            ddy_hjthkt = hk(1, jt, *y, xi, nd, nt)*hk(0, kt, *z, xi, nd, nt);
            ddz_hjthkt = hk(0, jt, *y, xi, nd, nt)*hk(1, kt, *z, xi, nd, nt);
            Fy += fy[jt+kt*ntp1]*hjthkt;
            Fz += fz[jt+kt*ntp1]*hjthkt;
            ddy_Fy += fy[jt+kt*ntp1]*ddy_hjthkt;
            ddy_Fz += fz[jt+kt*ntp1]*ddy_hjthkt;
            ddz_Fy += fy[jt+kt*ntp1]*ddz_hjthkt;
            ddz_Fz += fz[jt+kt*ntp1]*ddz_hjthkt;
         }
      }
      J = ddy_Fy*ddz_Fz-ddz_Fy*ddy_Fz;
      *y -= (Fy*ddz_Fz-Fz*ddz_Fy)/J;
      *z -= (Fz*ddy_Fy-Fy*ddy_Fz)/J;
      if(Fy*Fy+Fz*Fz<eps) break;
      if(counter>100){
         *y = LARGEVAL;
         *z = LARGEVAL;
         break;
      }
   }
   *y = fabs(*y);
   *z = fabs(*z);
   if(*y*(*y)>1.0 || *z*(*z)>1.0){
      *y = LARGEVAL;
      *z = LARGEVAL;
   }
   if((*y)*(*y)<(*z)*(*z)){
      tmp = *y;
      *y = *z;
      *z = tmp;
   }
}

void outputdata(double *intpltdfx, double *intpltdfr, double *intpltdft,
                                 int n, double *ye, double *ze, double ymax){
   int np1 = n+1;
   double y, z;
   int i, j, k;
   FILE *fpfx, *fpfr, *fpft, *fpeqp;
   char inputfilename1[] = "IntpltdFx.dat";
   char inputfilename2[] = "IntpltdFr.dat";
   char inputfilename3[] = "IntpltdFt.dat";
   char inputfilename4[] = "EquilibriumPoints.dat";
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
            fprintf(fpfx, MYFMT" "MYFMT" "MYFMT"\n",
                              y*ymax, z*ymax, intpltdfx[j+k*np1]);
            fprintf(fpfr, MYFMT" "MYFMT" "MYFMT"\n",
                              y*ymax, z*ymax, intpltdfr[j+k*np1]);
            fprintf(fpft, MYFMT" "MYFMT" "MYFMT"\n",
                              y*ymax, z*ymax, intpltdft[j+k*np1]);
         }
         fprintf(fpfx, "\n");
         fprintf(fpfr, "\n");
         fprintf(fpft, "\n");
      }
   fclose(fpfx);
   fclose(fpfr);
   fclose(fpft);
   fpeqp = fopen(inputfilename4, "w");
      fprintf(fpeqp, "# ye ze\n");
      for(i=0; i<8; i++){
         if(ye[i]==LARGEVAL) break;
         fprintf(fpeqp, MYFMT" "MYFMT"\n", ye[i]*ymax, ze[i]*ymax);
      }
   fclose(fpeqp);
}
