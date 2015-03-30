#include <math.h>
#include <stdio.h>
#include "mex.h"
#define NEIMAX 20
#include "opnml_mex5_allocs.c"

/************************************************************

  ####     ##     #####  ######  #    #    ##     #   #
 #    #   #  #      #    #       #    #   #  #     # #
 #       #    #     #    #####   #    #  #    #     #
 #  ###  ######     #    #       # ## #  ######     #
 #    #  #    #     #    #       ##  ##  #    #     #
  ####   #    #     #    ######  #    #  #    #     #

************************************************************/

/* PROTOTYPES */
void isort2(unsigned long n, int arr[], int brr[]);
void mexFunction(int            nlhs,
                 mxArray       *plhs[],
                 int            nrhs,
                 const mxArray *prhs[])

{
   /* all pointers start off as NULL */
   int **ele,**nei=NULL,*numnei=NULL,*itheta=NULL,*ilist=NULL;
   double *dnei, *dele, *x,*y,xc,yc,dx,dy;
   int maxnei=0;
   int i,j,k,kk,l,ll,ne,nn,no;

/* ---- check I/O arguments ----------------------------------------- */
   if (nrhs != 3)  mexErrMsgTxt("ele2nei requires 3 input arguments.");
   else if (nlhs!= 1)  mexErrMsgTxt("ele2nei requires 1 output arguments.");

/* ---- dereference input ------------------------------------------- */
   dele=mxGetPr(prhs[0]);
   ne=mxGetM(prhs[0]); 
   x=mxGetPr(prhs[1]);
   y=mxGetPr(prhs[2]);
   nn=mxGetM(prhs[1]);
   fprintf(stderr,"NE = %d\n",ne);  
   fprintf(stderr,"NN = %d\n",nn);  

/* ---- allocate space for int representation of dele &
        convert double element representation to int         -------- */
   ele=(int **)mxImatrix(1,ne,1,3);
   for (i=1;i<=ne;i++){
      for(j=1;j<=3;j++){
         ele[i][j]=((int)dele[(i-1)+(j-1)*ne]);
      }
   }
            
/* Allocate neighbor list array with matlab alloc routines in NRC style */
   nei=(int **)mxImatrix(1,nn,1,NEIMAX);
   numnei=(int *)mxIvector(1,nn);
     
   fprintf(stderr,"Computing neighbor list (Max=%d)...\n",NEIMAX);  
   for(k=1;k<=nn;k++){
      nei[k][1]=k;
      numnei[k]=1;
      for(l=1;l<=ne;l++){
         for(ll=1;ll<=3;ll++){
            if(ele[l][ll]==k)
               goto l83; 
         }                   
         goto l81;
 l83:    no=numnei[k];
         for(ll=1;ll<=3;ll++){
            for(kk=1;kk<=no;kk++)
               if(ele[l][ll]==nei[k][kk])goto l84;
            /* new neighbor */
            no++;
            nei[k][no]=ele[l][ll];
            numnei[k]=no;
 l84:       continue;
         }
 l81:    continue;
      }
   }

   maxnei=0;
   for(k=1;k<=nn;k++)
      maxnei=IMAX(maxnei,numnei[k]);
   maxnei--;
   fprintf(stderr,"Max Number of neighbors = %d\n",maxnei);

/* sort each node's neighbors into CCW around center node; ---------- */
   puts("Sorting neighboring node order on ANGLE");
   itheta=(int *)mxIvector(1,NEIMAX);
   ilist=(int *)mxIvector(1,NEIMAX);
   
   for(k=1;k<=nn;k++){
      xc=x[k-1];yc=y[k-1];
      for(j=2;j<=numnei[k];j++){         
         dx=x[nei[k][j]-1]-xc;
         dy=y[nei[k][j]-1]-yc;
         itheta[j-1]=(int)floor((atan2(dy,dx))*180./3.14159);
         if(itheta[j-1]<0)itheta[j-1]=itheta[j-1]+360;
         ilist[j-1]=nei[k][j];
      } 

      isort2(numnei[k]-1,itheta,ilist);
      for(j=2;j<=numnei[k];j++)
         nei[k][j]=ilist[j-1]; 

   } 
   puts("Sorting finished");
   
/* allocate and fill double neighbor vector for return */
   dnei=(double *)mxDvector(0,nn*maxnei);   
   for(j=0;j<maxnei;j++)
      for(k=0;k<nn;k++)
         dnei[k+j*nn]=(double)nei[k+1][j+2];
   
   plhs[0]=mxCreateDoubleMatrix(nn,maxnei,mxREAL); 
   mxSetPr(plhs[0],dnei);
          
/* ---- No need to free memory allocated with "mxCalloc"; MATLAB 
   does this automatically.  The CMEX allocation functions in 
   "opnml_allocs.c" use mxCalloc. ----------------------------------- */ 
   return;    


}
/*  NRC isort2 routine */

#define M 7
#define NSTACK 50
#define SWAP(AA,BB) temp=(AA);(AA)=(BB);(BB)=temp;

void isort2(unsigned long n, int arr[], int brr[])
{
   unsigned long i,ir=n,j,k,l=1;
   int *istack,jstack=0;
   int a,b,temp;

   istack=(int *)Ivector(1,NSTACK);
   for (;;) {
      if (ir-l < M) {
         for (j=l+1;j<=ir;j++) {
            a=arr[j];
            b=brr[j];
            for (i=j-1;i>=1;i--) {
               if (arr[i] <= a) break;
               arr[i+1]=arr[i];
               brr[i+1]=brr[i];
            }
            arr[i+1]=a;
            brr[i+1]=b;
         }
         if (!jstack) {
            free_Ivector(istack,1,NSTACK);
            return;
         }
         ir=istack[jstack];
         l=istack[jstack-1];
         jstack -= 2;
      } else {
         k=(l+ir) >> 1;
         SWAP(arr[k],arr[l+1])
         SWAP(brr[k],brr[l+1])
         if (arr[l+1] > arr[ir]) {
            SWAP(arr[l+1],arr[ir])
            SWAP(brr[l+1],brr[ir])
         }
         if (arr[l] > arr[ir]) {
            SWAP(arr[l],arr[ir])
            SWAP(brr[l],brr[ir])
         }
         if (arr[l+1] > arr[l]) {
            SWAP(arr[l+1],arr[l])
            SWAP(brr[l+1],brr[l])
         }
         i=l+1;
         j=ir;
         a=arr[l];
         b=brr[l];
         for (;;) {
            do i++; while (arr[i] < a);
            do j--; while (arr[j] > a);
            if (j < i) break;
            SWAP(arr[i],arr[j])
            SWAP(brr[i],brr[j])
         }
         arr[l]=arr[j];
         arr[j]=a;
         brr[l]=brr[j];
         brr[j]=b;
         jstack += 2;
         if (jstack > NSTACK) nrerror("NSTACK too small in isort2.");
         if (ir-i+1 >= j-l) {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
         } else {
            istack[jstack]=j-1;
            istack[jstack-1]=l;
            l=i;
         }
      }
   }
}
#undef M
#undef NSTACK



