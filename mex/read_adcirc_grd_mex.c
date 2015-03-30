/*

read_adcirc_grd_mex.c
Mex read function for ADCIRC grid files

the calling function (read_adcirc_grd.m) handles most error checking.  

7 Feb, 2011

*/

#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "opnml_mex5_allocs.c"
#define BUFFER_SIZE 132
#define MOD(i,j) (double)i/((double) Stride) - i/Stride
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

/* PROTOTYPES */
/* none */

/************************************************************

  ####     ##     #####  ######  #    #    ##     #   #
 #    #   #  #      #    #       #    #   #  #     # #
 #       #    #     #    #####   #    #  #    #     #
 #  ###  ######     #    #       # ## #  ######     #
 #    #  #    #     #    #       ##  ##  #    #     #
  ####   #    #     #    ######  #    #  #    #     #

************************************************************/

/* ---- read_adcirc_grd_mex will be called from :
        D=read_adcirc_grd(fname,debug); 
	where fname is the char string containing the grd filename
	---------------------------- */

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
   const char **fieldnames;
   int *dims;
   
   int i,j,it,k,junk;
   double *NND, *NNE;
   int NN, NE;

   int nfil;
   double defval;

   char *fortfilename,*firstline;
   int strlen,status;
   double *fortunit=NULL;
   FILE *fortfid=NULL, *fopen();
   double NaN=mxGetNaN();
   double FillThresh=-9999.;
   int SetFill=1;

   mxArray *nn, *ne;
   double *x, *y, *z, *e;
   mxArray *xx, *yy, *zz, *ee;

   /* --- check I/O arguments ------------------------------------------------------------ */
   if      (nrhs != 2)  mexErrMsgTxt("read_adcirc_grd_mex requires 2 input arguments.");
   else if (nlhs != 1)  mexErrMsgTxt("read_adcirc_grd_mex requires 1 output argument.");
      
   /* --- First input must be a string -------------------------------------------------- */
   if (mxIsChar(prhs[0]) != 1)mexErrMsgTxt("First input must be a string.");
   strlen=mxGetN(prhs[0])+1;
   fortfilename=mxCalloc(strlen,sizeof(char));
   status=mxGetString(prhs[0],fortfilename,strlen);
   if (status!=0)fprintf(stderr,"Input filename string extraction failed in READ_ADCIRC_FORT_MEX.");

   /* --- get debug flag ----------------------------------------------------------------- */
   double *temp;
   int debug;
   temp=mxGetPr(prhs[1]);
   debug=(int)temp[0];

   if (debug>1){
      printf("\nInside Mex file...\n");
      printf("\nnrhs          = %d\n",nrhs);
      printf("Debug         = %d\n",debug);
      printf("fortfilename=%s\n",fortfilename);
   }

   int nfields;
   nfields=8;
   fieldnames=mxCalloc(nfields,sizeof(*fieldnames));
   fieldnames[0]="filename";
   fieldnames[1]="name";
   fieldnames[2]="nn";
   fieldnames[3]="ne";
   fieldnames[4]="x";
   fieldnames[5]="y";
   fieldnames[6]="z";
   fieldnames[7]="e";

   if (!(fortfid = fopen(fortfilename,"r"))){
      printf("Open of %s failed.  Check for file. Terminal.\n",fortfilename);
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL); 
      return;
   }

   /* --- create return structure ----------------------------- */
   plhs[0] = mxCreateStructMatrix(1, 1, nfields , fieldnames);
   
   mxSetFieldByNumber(plhs[0], 0, 0, mxCreateString(fortfilename));
 
   firstline=mxCalloc(BUFFER_SIZE,sizeof(char));
   fgets(firstline,BUFFER_SIZE,fortfid);  
   if (debug>1)printf("\nfort.14 first line=%s\n",firstline);
   mxSetFieldByNumber(plhs[0], 0, 1, mxCreateString(firstline));

   NND=(double *) mxDvector(0,1);
   NNE=(double *) mxDvector(0,1);
   fscanf(fortfid,"%lf %lf",NNE,NND);
   NN=(int)NND[0];
   NE=(int)NNE[0];
   if(debug>1){
      printf("Number of ADCIRC nodes                 = %d\n",  NN); 
      printf("Number of ADCIRC elements              = %d\n",  NE); 
   }

   nn=mxCreateDoubleMatrix(1,1,mxREAL);
   mxFree(mxGetPr(nn));
   mxSetPr(nn,NND);
   mxSetFieldByNumber(plhs[0], 0, 2, nn);   

   ne=mxCreateDoubleMatrix(1,1,mxREAL);
   mxFree(mxGetPr(ne));
   mxSetPr(ne,NNE);
   mxSetFieldByNumber(plhs[0], 0, 3, ne);   

   x=   (double *) mxDvector(0,NN);
   y=   (double *) mxDvector(0,NN);
   z=   (double *) mxDvector(0,NN);
   if (debug>0)printf("Scanning %d Nodes ...\n",NN);
   for (i=0;i<NN;i++){
       fscanf(fortfid,"%d %lf %lf %lf",&junk,&x[i],&y[i],&z[i]);
   }

   if(debug>1)printf("Writing x structure ... \n");	   
   xx=mxCreateDoubleMatrix(NN,1,mxREAL);
   mxFree(mxGetPr(xx));
   mxSetPr(xx,x);
   mxSetFieldByNumber(plhs[0], 0, 4, xx);
   if(debug>1)printf("Writing y structure ... \n");	   
   yy=mxCreateDoubleMatrix(NN,1,mxREAL);
   mxFree(mxGetPr(yy));
   mxSetPr(yy,y);
   mxSetFieldByNumber(plhs[0], 0, 5, yy);
   if(debug>1)printf("Writing z structure ... \n");	   
   zz=mxCreateDoubleMatrix(NN,1,mxREAL);
   mxFree(mxGetPr(zz));
   mxSetPr(zz,z);
   mxSetFieldByNumber(plhs[0], 0, 6, zz);

   e = (double *) mxDvector(0,3*NE-1);
   if (debug>0)printf("Scanning %d Elements ...\n",NE);
   for (i=0;i<NE;i++){
       int k1,k2,k3,e1,e2,e3;
       k1=i+0*NE;
       k2=i+1*NE;
       k3=i+2*NE;
       fscanf(fortfid,"%d %d",&junk,&junk);
       fscanf(fortfid,"%d %d %d",&e1,&e2,&e3);
       e[k1]=e1; e[k2]=e2; e[k3]=e3;
   }
   if(debug>1)printf("Writing e structure ... \n");	   
   ee=mxCreateDoubleMatrix(NE,3,mxREAL);
   mxFree(mxGetPr(ee));
   mxSetPr(ee,e);
   mxSetFieldByNumber(plhs[0], 0, 7, ee);




   fclose(fortfid);
   if (debug>1){
      printf("\nLeaving Mex file...\n");
   }

   return;
}
