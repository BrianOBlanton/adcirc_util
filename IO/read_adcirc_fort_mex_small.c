/*

read_adcirc_fort_mex_small.c
Direct read function for ADCIRC fort files 

the calling function (read_adcirc_fort.m) handles
most error checking.  

19 Dec, 2007

Reduced functionality for debugging

*/

#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "opnml_mex5_allocs.c"
#define BUFFER_SIZE 132

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


void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
   
   const char **fieldnames;
   int i,j,it,k; 
   double *time;
   char *filename,*header,*linejunk;
   int strlen,status;
   FILE *fid=NULL, *fopen();
   double NaN=mxGetNaN();
   double FillThresh=-9999.;
   int NDSETS,NSTEP,IFLAG,INN;
   double *DT,*NN,*z;

   int debug=0;
   int nfields=6;
   
   mxArray *dt, *nn, *ttime, *zz;
      
   /* --- check I/O arguments ------------------------------------------------------------ */
   if      (nrhs != 1)  mexErrMsgTxt("read_adcirc_fort_mex_small requires 1 input arguments.");
   else if (nlhs != 1)  mexErrMsgTxt("read_adcirc_fort_mex_small requires 1 output argument.");
      
   /* --- Input must be a string --------------------------------------------------------- */
   if (mxIsChar(prhs[0]) != 1)mexErrMsgTxt("First input must be a string.");

   if (debug>0){
      fprintf(stdout,"\nDebug       = %d\n",debug);
      fprintf(stdout,"\nnrhs        = %d\n",nrhs);
   }

   /* --- open data file ---------------------------------------------------------------- */
   strlen=mxGetN(prhs[0])+1;
   filename=mxCalloc(strlen,sizeof(char));
   status=mxGetString(prhs[0],filename,strlen);
   if (status!=0)fprintf(stderr,"Input filename string extraction failed in READ_ADCIRC_FORT_MEX.");

   if (debug>0)printf("\nfilename=%s\n",filename);
 
   // if open fails, return empty
   if (!(fid = fopen(filename,"r"))){
      printf("Open of %s failed.  Check for file. Terminal.\n",filename);
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL); 
      return;}
   
   // set return struct field names
   fieldnames=mxCalloc(nfields,sizeof(*fieldnames));
   fieldnames[0]="filename";
   fieldnames[1]="header";
   fieldnames[2]="dt";
   fieldnames[3]="nn";
   fieldnames[4]="t";
   fieldnames[5]="zeta";

   /* --- create return structure ----------------------------- */
   plhs[0] = mxCreateStructMatrix(1, 1, nfields , fieldnames);
 
   header=mxCalloc(BUFFER_SIZE,sizeof(char));
   linejunk=mxCalloc(BUFFER_SIZE,sizeof(char));
   fgets(header,BUFFER_SIZE,fid);  
   if (debug>0)printf("\nfort.XX header line=%s\n",header);

   DT=(double *) mxDvector(0,1);
   NN=(double *) mxDvector(0,1);

   fscanf(fid,"%d %lf %lf %d %d",&NDSETS,NN,DT,&NSTEP,&IFLAG);
   fgets(linejunk,BUFFER_SIZE,fid);  
   INN=(int)NN[0];

   if (debug>0){
      printf("Expected Number of Data Sets in file   = %d\n",  NDSETS); 
      printf("Number of nodes                        = %d\n",  INN); 
      printf("Time Step in fort file                 = %lf\n", DT[0]); 
      printf("Number of Model steps per output       = %d\n",  NSTEP);
   }

   /* --- (Partly) Fill return structure ---------------------- */
   dt=mxCreateDoubleMatrix(1,1,mxREAL);
   mxFree(mxGetPr(dt));
   mxSetPr(dt,DT);
   mxSetFieldByNumber(plhs[0], 0, 2, dt);

   nn=mxCreateDoubleMatrix(1,1,mxREAL);
   mxFree(mxGetPr(nn));
   mxSetPr(nn,NN);
	   
   mxSetFieldByNumber(plhs[0], 0, 0, mxCreateString(filename));
   mxSetFieldByNumber(plhs[0], 0, 1, mxCreateString(header));
   mxSetFieldByNumber(plhs[0], 0, 2, dt);
   mxSetFieldByNumber(plhs[0], 0, 3, nn);   
   
   /* --- Scan data part of fort.?? file ---------------------- */
   z=   (double *) mxDvector(0,INN*NDSETS);
   time=(double *) mxDvector(0,NDSETS);

   /* init to NaN */
   for (j=0;j<INN*NDSETS;j++){z[j]=NaN;} 
   for (j=0;j<NDSETS;j++){time[j]=NaN;}  

   for (i=0;i<NDSETS;i++){
      fscanf(fid,"%lf %d",&time[i],&it);
      if (debug>1)printf("Reading!! %d %lf %d\n",i,time[i],it); 
      for (j=0;j<INN;j++){fscanf(fid,"%d %lf",&k,&z[j+i*INN]);}
      for (j=0;j<INN;j++){
         if (z[j+i*INN]<FillThresh){z[j+i*INN]=NaN;} ;
      }   
   }

   // fill data parts of return struct
   
   if(debug>0)printf("\nWriting time to structure ... \n");	   
   ttime=mxCreateDoubleMatrix(1,NDSETS,mxREAL);
   mxFree(mxGetPr(ttime));
   mxSetPr(ttime,time);
   mxSetFieldByNumber(plhs[0], 0, 4, ttime);

   if(debug)fprintf(stdout,"\nWriting data to structure %dx%d ... \n",INN,NDSETS);	   
   zz=mxCreateDoubleMatrix(INN,NDSETS,mxREAL);
   mxFree(mxGetPr(zz));
   mxSetPr(zz,z);
   mxSetFieldByNumber(plhs[0], 0, 5, zz);

   // Destroy/free temp arrays
   
   fclose(fid);
   return;
}
  
