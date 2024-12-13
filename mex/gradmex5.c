#include <math.h> 
#include <stdio.h> 
#include "mex.h" 
#include "opnml_mex5_allocs.c" 



/* PROTOTYPES */
#ifdef __STDC__
void bandmsolve_c(int,
                  double **,
                  double *,
                  int,
                  int);
void grad(int nn,
          int ne,
          double *x,
          double *y,
          int *ele,
          double *q,
          double *grax,
          double *gray);
#endif

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
   double *x, *y, *q;
   double *dele;
   double *grax, *gray;

/* ---- gradmex will be called as :
        [gradx,grady]=gradmex(x,y,ele,q); 
                                      ------------------------------- */
        
/* ---- check I/O arguments ----------------------------------------- */
   if (nrhs != 4) 
      mexErrMsgTxt("gradmex requires 4 input arguments.");
   else if (nlhs != 2) 
      mexErrMsgTxt("gradmex requires 2 output arguments.");

/* ---- dereference input arrays ------------------------------------ */
   x=mxGetPr(prhs[0]);
   y=mxGetPr(prhs[1]);
   dele=mxGetPr(prhs[2]);
   q=mxGetPr(prhs[3]);
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
   
/* ---- allocate space for gradient lists  following 
        NRC allocation style                            ------------- */
   grax = (double *)  mxDvector(0,nn-1);
   gray = (double *)  mxDvector(0,nn-1);
   
/* ---- call gradient routine  -------------------------------------- */
   grad(nn,ne,x,y,ele,q,grax,gray);
     
/* ---- Set elements of return matricies, 
        pointed to by plhs[0] & plhs[1] 
                                        ---------------------------- */
   plhs[0]=mxCreateDoubleMatrix(nn,1,mxREAL); 
   mxSetPr(plhs[0],grax);
   plhs[1]=mxCreateDoubleMatrix(nn,1,mxREAL); 
   mxSetPr(plhs[1],gray);

/* --- No need to free memory allocated with "mxCalloc"; MATLAB 
   does this automatically.  The CMEX allocation functions in 
   "opnml_allocs.c" use mxCalloc. ----------------------------------- */    
    
   return;
}

/*----------------------------------------------------------------------

  ####   #####     ##    #####
 #    #  #    #   #  #   #    #
 #       #    #  #    #  #    #
 #  ###  #####   ######  #    #
 #    #  #   #   #    #  #    #
  ####   #    #  #    #  #####
 
----------------------------------------------------------------------*/

void grad(int nn,int ne,
          double *x, 
          double *y,
          int *ele,
          double *q, 
          double *grax,
          double *gray)
#define ELE(i,j,m) ele[i+m*j]
#define DZILCH (double)0.

{
   double *ar,**av,**dx,**dy;
   double area6,area12;
   int i,j,l,m,n;
   int nbw=0, nh=0, nhm=0, l01, l02, l12;

/* ---- Compute half bandwidth -------------------------------------- */
   for (i=0;i<ne;i++){
      l01=abs(ELE(i,0,ne)-ELE(i,1,ne));
      l02=abs(ELE(i,0,ne)-ELE(i,2,ne));
      l12=abs(ELE(i,1,ne)-ELE(i,2,ne));
      nhm=IMAX(l01,IMAX(l02,l12));  
      nh=IMAX(nh,nhm);    
   }
   nbw=2*nh+1;
   //printf("BW = %d\n",nbw);

/* ---- ALLOCATE space for dx, dy, ar, av and dv -------------------- */
   ar = (double *)  mxDvector(0,ne-1);
   av = (double **) mxDmatrix(0,nbw-1,0,nn-1);
   dx = (double **) mxDmatrix(0,2,0,ne-1);
   dy = (double **) mxDmatrix(0,2,0,ne-1);   
   //printf("Local Arrays Allocated\n");
            
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
   //printf("Element areas and edge lengths computed\n");
   
/* ---- Assemble LHS (av) ------------------------------------------- */
   for(l=0;l<ne;l++){
      area6=ar[l]/(double)6.;
      area12=ar[l]/(double)12.;
      for(j=0;j<3;j++){
         m=ELE(l,j,ne);
         for(i=0;i<3;i++){
            n=nh+ELE(l,i,ne)-ELE(l,j,ne);
            if(i==j) av[n][m]+=area6;
            else     av[n][m]+=area12;      
         }
      }
   }
   //puts("LHS (av) assembled");

/* ---- Assemble RHS (grax,gray) ------------------------------------ */
   for(l=0;l<ne;l++){
      for(j=0;j<3;j++){
         m=ELE(l,j,ne);
         for(i=0;i<3;i++){
            grax[m]+=q[ELE(l,i,ne)]*dy[i][l]/(double)6.;
            gray[m]-=q[ELE(l,i,ne)]*dx[i][l]/(double)6.;
         }
      }
   }
   //puts("RHS (grax,gray) assembled");

/* ---- Call banded Matrix solver to triangularize av --------------- */
   bandmsolve_c(0,av,grax,nn,nh);
   //puts("Left side triangularized");

/* ---- Call banded Matrix solver to solve for (grax,gray) ---------- */
   bandmsolve_c(1,av,grax,nn,nh);
   bandmsolve_c(1,av,gray,nn,nh);
   //puts("Gradient computed");
 
   return;     
}



