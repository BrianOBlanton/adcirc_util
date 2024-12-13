#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "opnml_mex5_allocs.c"

/* PROTOTYPES */
void bandmsolve_c(int,
                  double **,
                  double *,
                  int,
                  int);
void curl(int nn,
          int ne,
          double *x,
          double *y,
          int *ele,
          double *u,
          double *v,
          double *dv);

/************************************************************

  ####     ##     #####  ######  #    #    ##     #   #
 #    #   #  #      #    #       #    #   #  #     # #
 #       #    #     #    #####   #    #  #    #     #
 #  ###  ######     #    #       # ## #  ######     #
 #    #  #    #     #    #       ##  ##  #    #     #
  ####   #    #     #    ######  #    #  #    #     #

************************************************************/


   void mexFunction(int            nlhs,
                    mxArray       *plhs[],
                    int            nrhs,
                    const mxArray *prhs[])
{
   int cnt,*ele,i,j,nn,ne;
   double *x, *y, *u, *v;
   double *dele;
   double *vt;

/* ---- curlmex will be called as :
        vt=curlmex(x,y,ele,u,v); 
                                      ------------------------------- */
        
/* ---- check I/O arguments ----------------------------------------- */
   if (nrhs != 5) 
      mexErrMsgTxt("curlmex5 requires 5 input arguments.");
   else if (nlhs != 1) 
      mexErrMsgTxt("curlmex5 requires 1 output arguments.");

/* ---- dereference input arrays ------------------------------------ */
   x=mxGetPr(prhs[0]);
   y=mxGetPr(prhs[1]);
   dele=mxGetPr(prhs[2]);
   u=mxGetPr(prhs[3]);
   v=mxGetPr(prhs[4]);
   nn=mxGetM(prhs[0]);
   ne=mxGetM(prhs[2]);   

/* ---- allocate space for int representation of dele &
        convert double element representation to int  &
        shift node numbers toward 0 by 1 for proper indexing -------- */
   ele=(int *)mxIvector(0,3*ne);
   for (i=0;i<3*ne;i++){
      ele[i]=(int)dele[i];
      ele[i]=ele[i]-1;
   }
   
/* ---- allocate space for divergence list dv following 
        NRC allocation style                            ------------- */
   vt = (double *)  mxDvector(0,nn-1);
   
/* ---- call curl routine  ------------------------------------ */
   curl(nn,ne,x,y,ele,u,v,vt);
     
/* ---- Set elements of return matrix, pointed to by plhs[0] -------- */
   plhs[0]=mxCreateDoubleMatrix(nn,1,mxREAL); 
   mxSetPr(plhs[0],vt);   

/* --- No need to free memory allocated with "mxCalloc"; MATLAB 
   does this automatically.  The CMEX allocation functions in 
   "opnml_allocs.c" use mxCalloc. ----------------------------------- */    
    
   return;
}

/*----------------------------------------------------------------------
 
----------------------------------------------------------------------*/
void curl(int nn,int ne,
             double *x, 
             double *y,
             int *ele,
             double *u, 
             double *v,
             double *vt)
#define ELE(i,j,m) ele[i+m*j]

{
   double *ar,**av,**dx,**dy;
   double area6,area12;
   int i,j,l,m,n;
   int nbw=0, nh=0, nhm=0, n0 ,n1, n2, l01, l02, l12;

/* ---- Compute half bandwidth -------------------------------------- */
   for (i=0;i<ne;i++){
      l01=abs(ELE(i,0,ne)-ELE(i,1,ne));
      l02=abs(ELE(i,0,ne)-ELE(i,2,ne));
      l12=abs(ELE(i,1,ne)-ELE(i,2,ne));
      nhm=IMAX(l01,IMAX(l02,l12));  
      nh=IMAX(nh,nhm);    
   }
   nbw=2*nh+1;
   fprintf(stderr,"BW = %d\n",nbw);

/* ---- ALLOCATE space for dx, dy, ar, av and dv -------------------- */
   ar = (double *)  mxDvector(0,ne-1);
   av = (double **) mxDmatrix(0,nbw-1,0,nn-1);
   dx = (double **) mxDmatrix(0,2,0,ne-1);
   dy = (double **) mxDmatrix(0,2,0,ne-1);   
   puts("Local Arrays Allocated");
            
/* ---- Compute dx,dy and element areas ----------------------------- */
   for(l=0;l<ne;l++){
      dx[0][l]=x[ELE(l,1,ne)]-x[ELE(l,2,ne)];
      dx[1][l]=x[ELE(l,2,ne)]-x[ELE(l,0,ne)];
      dx[2][l]=x[ELE(l,0,ne)]-x[ELE(l,1,ne)];
      dy[0][l]=y[ELE(l,1,ne)]-y[ELE(l,2,ne)];
      dy[1][l]=y[ELE(l,2,ne)]-y[ELE(l,0,ne)];
      dy[2][l]=y[ELE(l,0,ne)]-y[ELE(l,1,ne)];
      ar[l]=.5*(x[ELE(l,0,ne)]*dy[0][l]+
                x[ELE(l,1,ne)]*dy[1][l]+
                x[ELE(l,2,ne)]*dy[2][l]);      
   }
   puts("Element areas and edge lengths computed");
   
/* ---- Assemble LHS (av) ------------------------------------------- */
   for(l=0;l<ne;l++){
      area6=ar[l]/6.;
      area12=ar[l]/12.;
      for(j=0;j<3;j++){
         m=ELE(l,j,ne);
         for(i=0;i<3;i++){
            n=nh+ELE(l,i,ne)-ELE(l,j,ne);
            if(i==j) av[n][m]+=area6;
            else     av[n][m]+=area12;      
         }
      }
   }
   puts("LHS (av) assembled");

/* ---- Assemble RHS (vt) ------------------------------------------- */
   for(l=0;l<ne;l++){
      for(j=0;j<3;j++){
         m=ELE(l,j,ne);
         for(i=0;i<3;i++)
              vt[m]+=v[ELE(l,i,ne)]*dy[i][l]/6.+u[ELE(l,i,ne)]*dx[i][l]/6.;
      }
   }
   puts("RHS (vt) assembled");

/* ---- Call banded Matrix solver to triangularize av --------------- */
   bandmsolve_c(0,av,vt,nn,nh);
   puts("Left side triangularized");
   
/* ---- Call banded Matrix solver to solve for dv ------------------- */
   bandmsolve_c(1,av,vt,nn,nh);
   puts("Curl computed");
 
   return;     
}

