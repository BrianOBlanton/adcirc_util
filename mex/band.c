#include <stdio.h>
#include "opnml_mex5_allocs.c"
/*********  BANDMSOLVE SUBROUTINE   ******************/

void bandmsolve_c(int kkk, double **b, double *r, int neq, int hbw)
#define DZILCH (double)0.
{
/* Asymmetric band matrix equation solver
   doctored to ignore 0's in the lu decomp step
   
   Input kkk = 0 ==>>  triangularize the band matrix b
         kkk = 1 ==>>  solve for right side r, return solution in r
*/
   
   int bw,i,ib,ihbp,j,jc,ji,jp,k;
   int kc,ki,kk,kr,lim,mr,nn;
   double c,sum,pivot;
   
/* Compute BW */
   bw=2*hbw+1;

   ihbp=hbw+1;  
   
   if (kkk==0){

/* ---- Triangularize b using Doolittle Method ---- */
      for (k=1;k<=neq-1;++k) {
         pivot=b[ihbp-1][k-1];
         kk=k+1;
         kc=ihbp;
         for (i=kk;i<=neq;++i){
            --kc;
            if(kc<=0)goto L10;
            c=-b[kc-1][i-1]/pivot;
            if (c==(float)0)goto L21;
            b[kc-1][i-1]=c;
            ki=kc+1;
            lim=kc+hbw;
            for (j=ki;j<=lim;++j){
               jc=ihbp+j-kc;
               b[j-1][i-1]+=c*b[jc-1][k-1];
            }
L21:     ;
         }
L10:    ;}
   }
   
   else if (kkk==1){

/* ---- Modify load vector ---- */
      nn=neq+1;
      for (i=2;i<=neq;++i){
         jc=ihbp-i+1;
         ji=1;
         if (jc<=0){
            jc=1;
            ji=i-ihbp+1;
         }
         sum=DZILCH;
         for (j=jc;j<=hbw;++j){
            sum+=b[j-1][i-1]*r[ji-1];
            ++ji;
         }
         r[i-1]+=sum;
      }

/* ---- Back Solution ---- */
      r[neq-1]/=b[ihbp-1][neq-1];
      for (ib=2;ib<=neq;++ib){
         i=nn-ib;
         jp=i;
         kr=ihbp+1;
         mr=IMIN(bw,hbw+ib);
         sum=DZILCH;
         for (j=kr;j<=mr;++j){
            ++jp;
            sum+=b[j-1][i-1]*r[jp-1];
         }
         r[i-1]=(r[i-1]-sum)/b[ihbp-1][i-1];
      }
   }
   else{
      fprintf(stderr,"Incorrect flag passed in 1st argument to bandmsolve.\n");
      exit(0);
   }
}
