/*

read_adcirc_fort_compact_mex.c
Mex read function for ADCIRC fort files, for the "compacted" global output format

the calling function (read_adcirc_fort.m) handles
most error checking.  

12 Jul, 2007

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

/* ---- read_adcirc_fort_compact_mex will be called from :
        D=read_adcirc_fort(fname,unit,debug,stride); 
	where fname is the char string containing the fort.?? filename
	and unit is the ?? number.  This selects what gets read.
	---------------------------- */

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
   const char **fieldnames;
   int *dims;
   
   int i,j,it,k; 
   double *time, *iter;

   int nfil;
   double defval;

   char *fortfilename,*header,*linejunk;
   int strlen,status,*supported_units=NULL,nunits,ifortunit;
   double *fortunit=NULL;
   FILE *fortfid=NULL, *fopen();
   double NaN=mxGetNaN();
   double FillThresh=-9999.;
   int SetFill=1;
   int itest;
         
   int NDSETS1,NDSETS2,NSTEP,IRTYPE,NN,ndsets;
   int actual_steps,setcnt;
   double *DT,*NND;
   mxArray *dt, *nn, *ttime, *iiter;
   mxArray *zz,*uu,*vv;
   double *z,ztempsingle,*u,*v,utempsingle,vtempsingle;
   double stridetest;
   double timestrip;
   double *temp;
   int debug;
   int Stride;
   int Strip;
   int IterStart;
   int IterEnd;
   int ii;
   int nfields;
     
   nunits=9;
   supported_units=mxIvector(0,nunits-1);
   supported_units[0]=61;
   supported_units[1]=62;
   supported_units[2]=63;
   supported_units[3]=64;
   supported_units[4]=69;
   supported_units[5]=71;
   supported_units[6]=72;
   supported_units[7]=73;
   supported_units[8]=74;

   /* --- check I/O arguments ------------------------------------------------------------ */
   if      (nrhs != 7)  mexErrMsgTxt("read_adcirc_fort_mex_compact requires 7 input arguments.");
   else if (nlhs != 1)  mexErrMsgTxt("read_adcirc_fort_mex_compact requires 1 output argument.");
      
   /* --- First input must be a string -------------------------------------------------- */
   if (mxIsChar(prhs[0]) != 1)mexErrMsgTxt("First input must be a string.");

   /* --- Second input must be from set={63,64,...}  ------------------------------------- */
   fortunit=mxGetPr(prhs[1]);
   ifortunit=(int)fortunit[0];
   
   /* --- get debug flag ----------------------------------------------------------------- */
   temp=mxGetPr(prhs[2]);
   debug=(int)temp[0];

   /* --- get stride flag ---------------------------------------------------------------- */
   temp=mxGetPr(prhs[3]);
   Stride=(int)temp[0];

   /* --- get strip flag ---------------------------------------------------------------- */
   temp=mxGetPr(prhs[4]);
   Strip=(int)temp[0];

   /* --- get starting iter number ------------------------------------------------------ */
   temp=mxGetPr(prhs[5]);
   IterStart=(int)temp[0];
   
   /* --- get ending iter number -------------------------------------------------------- */
   temp=mxGetPr(prhs[6]);
   IterEnd=(int)temp[0];
   
   if (debug>0){
      printf("\nInside Mex file...\n");
      printf("\nnrhs          = %d\n",nrhs);
      printf("Debug         = %d\n",debug);
      printf("Unit number   = %d\n",ifortunit);
      printf("Stride        = %d\n",Stride);
      printf("IterStart     = %d\n",IterStart);
      printf("IterEnd       = %d\n",IterEnd);      
      printf("Strip         = %d\n",Strip);
      printf("SetFill       = %d\n",SetFill);
   }

   /* check unit number */
   for (ii=0;ii<nunits;ii++){if (ifortunit == supported_units[ii])goto l20;}
   mexErrMsgTxt("Unsupported ADCIRC output unit number.");
   l20: if (debug>0) printf("unit index=%d\n",ii);

   /* --- open fort.?? file ------------------------------------------ */
   strlen=mxGetN(prhs[0])+1;
   fortfilename=mxCalloc(strlen,sizeof(char));
   status=mxGetString(prhs[0],fortfilename,strlen);
   if (status!=0)fprintf(stderr,"Input filename string extraction failed in READ_ADCIRC_FORT_MEX.");

   if (debug>0)printf("\nfortfilename=%s\n",fortfilename);
 
   if (!(fortfid = fopen(fortfilename,"r"))){
      printf("Open of %s failed.  Check for file. Terminal.\n",fortfilename);
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL); 
      return;
   }
   
   /* --- specify structure field names, ifortunit based -------------- */
   if ( (ifortunit == 61) | (ifortunit == 63) | (ifortunit == 69) ){      
      nfields=7;
      fieldnames=mxCalloc(nfields,sizeof(*fieldnames));
      fieldnames[0]="filename";
      fieldnames[1]="header";
      fieldnames[2]="dt";
      fieldnames[3]="nn";
      fieldnames[4]="t";
      fieldnames[5]="iter";
      fieldnames[6]="zeta";
   }
   else if ( (ifortunit == 71) | (ifortunit == 73) ){      
      nfields=7;
      fieldnames=mxCalloc(nfields,sizeof(*fieldnames));
      fieldnames[0]="filename";
      fieldnames[1]="header";
      fieldnames[2]="dt";
      fieldnames[3]="nn";
      fieldnames[4]="t";
      fieldnames[5]="iter";
      fieldnames[6]="pres";
   }
   else if( (ifortunit == 62) | (ifortunit == 64) ){
      nfields=8;
      fieldnames=mxCalloc(nfields,sizeof(*fieldnames));
      fieldnames[0]="filename";
      fieldnames[1]="header";
      fieldnames[2]="dt";
      fieldnames[3]="nn";
      fieldnames[4]="t";
      fieldnames[5]="iter";
      fieldnames[6]="ubar";
      fieldnames[7]="vbar";
   }
   else if( (ifortunit == 72) | (ifortunit == 74)){
      nfields=8;
      fieldnames=mxCalloc(nfields,sizeof(*fieldnames));
      fieldnames[0]="filename";
      fieldnames[1]="header";
      fieldnames[2]="dt";
      fieldnames[3]="nn";
      fieldnames[4]="t";
      fieldnames[5]="iter";
      fieldnames[6]="uwind";
      fieldnames[7]="vwind";
   }

   /* --- create return structure ----------------------------- */
   plhs[0] = mxCreateStructMatrix(1, 1, nfields , fieldnames);
 
   header=mxCalloc(BUFFER_SIZE,sizeof(char));
   linejunk=mxCalloc(BUFFER_SIZE,sizeof(char));
   fgets(header,BUFFER_SIZE,fortfid);  
   if (debug>0)printf("\nfort.XX header line=%s\n",header);

   DT=(double *) mxDvector(0,1);
   NND=(double *) mxDvector(0,1);

   fscanf(fortfid,"%d %lf %lf %d %d",&NDSETS1,NND,DT,&NSTEP,&IRTYPE);
   fgets(linejunk,BUFFER_SIZE,fortfid);  
  
   /* Set NDSETS2 to min of NDSETS1 and IterEnd */
   NDSETS2=IMIN(NDSETS1,IterEnd);
   
   ndsets=ceil((float)(NDSETS2-Strip)/((float)Stride));

   NN=(int)NND[0];

   if (debug>0){
      printf("Total Number of Data Sets in file      = %d\n",  NDSETS1); 
      printf("Max Number of Data Sets to read        = %d\n",  NDSETS2); 
      printf("Data Set Stride                        = %d\n",  Stride);
      printf("Data Set Strip                         = %d\n",  Strip);
      printf("Number of Data Sets to return          = %d\n",  ndsets);
      printf("Number of ADCIRC nodes                 = %d\n",  NN); 
      printf("Time Step in fort file                 = %lf\n", DT[0]); 
      printf("Number of Model steps per output       = %d\n",  NSTEP);
      if (SetFill){
         printf("Field with values less than %lf will be set to NaN\n",FillThresh);
      }
   }

   if (Strip>=NDSETS1){
      mexErrMsgTxt("Strip value cannot be larger than number of steps in input file. Lower your expectations.");
   }
   
   if (ndsets<1){
      mexErrMsgTxt("After stripping and striding, there is nothing left. Lower your expectations.");
   }
   
   /* --- (Partly) Fill return structure ---------------------- */
   dt=mxCreateDoubleMatrix(1,1,mxREAL);
   mxFree(mxGetPr(dt));
   mxSetPr(dt,DT);
   
   nn=mxCreateDoubleMatrix(1,1,mxREAL);
   mxFree(mxGetPr(nn));
   mxSetPr(nn,NND);
	   
   mxSetFieldByNumber(plhs[0], 0, 0, mxCreateString(fortfilename));
   mxSetFieldByNumber(plhs[0], 0, 1, mxCreateString(header));
   mxSetFieldByNumber(plhs[0], 0, 2, dt);
   mxSetFieldByNumber(plhs[0], 0, 3, nn);   
  
   /* --- switch on fort number --------------------------------------- */
   
   /* --- Case fort.61,63,73  ----------------------------------------- */
   if (ifortunit == 61 | ifortunit == 63 | ifortunit == 73 ){      	   

      /* --- allocate return arrays         ---------------------- */
      z=   (double *) mxDvector(0,NN*ndsets);
      time=(double *) mxDvector(0,ndsets);
      iter=(double *) mxDvector(0,ndsets);

      /* z is initialized below */
      for (j=0;j<ndsets;j++){time[j]=NaN;}  /* init to NaN */

      /* --- Scan data part of fort.?? file --------------------------- */
                 
      /* Strip off Strip timesteps */
      if (debug>0)printf("Stripping off %d time steps.\n",Strip);
      for (i=0;i<Strip;i++){
	      fscanf(fortfid,"%lf %d %d %lf",&timestrip,&it,&nfil,&defval);
	      for (j=0;j<nfil;j++){
	         fscanf(fortfid,"%d %lf",&k,&ztempsingle);
         }
         if (debug>0)printf("Skipped %d %lf %d %d\n",i,timestrip,it,nfil);
      } 
		        
      actual_steps=0;
      for (i=Strip;i<NDSETS2;i++){

         double firstval,thisval;
         
         itest=(i+1) % Stride;

         if ( itest == 0 ){
	    
            fscanf(fortfid,"%lf %lf %d %lf",&time[actual_steps],&iter[actual_steps],&nfil,&defval);
            
            if (debug>1)printf("Scanning %d %d %lf %d %d %lf \n",
                i,actual_steps,time[actual_steps],iter[actual_steps],nfil,defval); 

	    /* initialize z array to defval */
	    for (j=0;j<NN;j++){
	       int idx=j+actual_steps*NN;
               z[idx]=defval;
            } 
            
            for (j=0;j<nfil;j++){
	       fscanf(fortfid,"%d %lf",&k,&thisval);
	       z[k-1+actual_steps*NN]=thisval;
	       if (j == 0)firstval=thisval;
	    }
            
	    if (SetFill){
               for (j=0;j<NN;j++){
                  if (z[j+actual_steps*NN]<FillThresh){
                     z[j+actual_steps*NN]=NaN;
                  } 
	       }
	    }

            if (debug>2) printf("First val at this iteration :  %d : %lf\n",k,firstval);
            if (debug>2) printf("Last val at this iteration :  %d : %lf\n",k,thisval);
            if (debug>2) printf("Scanned %d values at it=%d\n",nfil,it);
	    	    
            actual_steps++;
         }

        /* Skip this iter */
         else{
            double timetemp, ztemp;
            fscanf(fortfid,"%lf %d %d %lf",&timetemp,&it,&nfil,&defval);
            if (debug>1) printf("Skipped !! %d %d\n",nfil,it);
            for (j=0;j<nfil;j++){
               fscanf(fortfid,"%d %lf",&k,&ztemp);
	    }
         }
      }

      if(debug>0)printf("Writing time structure ... \n");	   
      ttime=mxCreateDoubleMatrix(1,ndsets,mxREAL);
      mxFree(mxGetPr(ttime));
      mxSetPr(ttime,time);
      mxSetFieldByNumber(plhs[0], 0, 4, ttime);

      if(debug>0)printf("Writing iter structure ... \n");	   
      iiter=mxCreateDoubleMatrix(1,ndsets,mxREAL);
      mxFree(mxGetPr(iiter));
      mxSetPr(iiter,iter);
      mxSetFieldByNumber(plhs[0], 0, 5, iiter);
      
      if(debug>0)printf("Writing zeta structure %dx%d ... \n",NN,ndsets);	   
      zz=mxCreateDoubleMatrix(NN,ndsets,mxREAL);
      mxFree(mxGetPr(zz));
      mxSetPr(zz,z);
      mxSetFieldByNumber(plhs[0], 0, 6, zz);
   
   } 
   /* --- END CASE FORT.61,63,73 -------------------------------------- */

   /* --- Case fort.62,64,74 ------------------------------------------ */
	   
   else if(ifortunit == 62 | ifortunit == 64 | ifortunit == 72 | ifortunit == 74 ){
   
      /* --- allocate return arrays         ---------------------- */
      u=   (double *) mxDvector(0,NN*ndsets);
      v=   (double *) mxDvector(0,NN*ndsets);
      time=(double *) mxDvector(0,ndsets);
      iter=(double *) mxDvector(0,ndsets);

      for (j=0;j<NN*ndsets;j++){
         u[j]=NaN; /* init to NaN */
         v[j]=NaN;
      } 
      for (j=0;j<ndsets;j++){time[j]=NaN;}  /* init to NaN */
	   
      /* --- Scan data part of fort.?? file ---------------------- */

      /* Strip off Strip timesteps */
      if (debug>0)printf("Stripping off %d time steps.\n",Strip);
      for (i=0;i<Strip;i++){

	 fscanf(fortfid,"%lf %d %d %lf",&timestrip,&it,&nfil,&defval);
	 if (debug>0)printf("Skipping %lf %d %d\n",timestrip,it,nfil);
	 for (j=0;j<nfil;j++){
	   fscanf(fortfid,"%d %lf %lf",&k,&utempsingle,&vtempsingle);
         }
	 if (debug>1)printf("   Last val:  %lf %lf\n",utempsingle,vtempsingle);
      } 

      actual_steps=0;
      for (i=Strip;i<NDSETS2;i++){

         itest=(i+1) % Stride;

         if ( itest== 0 ){

            fscanf(fortfid,"%lf %lf %d %lf",&time[actual_steps],&iter[actual_steps],&nfil,&defval);
            
            if (debug>1)printf("Scanning  %d %d %lf %d %d %lf \n",
                  i,actual_steps,time[actual_steps],iter[actual_steps],nfil,defval); 
	    
	    /* initialize u,v arrays to defval */
            for (j=0;j<NN;j++){
	       int idx=j+actual_steps*NN;
	       u[idx]=defval;
 	       v[idx]=defval;
            }
           
            for (j=0;j<nfil;j++){
               fscanf(fortfid,"%d %lf %lf",&k,&utempsingle,&vtempsingle);
 	       u[k-1+actual_steps*NN]=utempsingle;
 	       v[k-1+actual_steps*NN]=vtempsingle;
            }
	    
            if (debug>2) printf("Last u,v at this iteration : %d : %lf,%lf\n",
                  k,utempsingle,vtempsingle);
            if (debug>2) printf("Scanned %d values at it=%d\n",nfil,it);

            actual_steps++;
         }

	 /* Skip this iter */
         else{
	    double timetemp, utemp, vtemp;
            fscanf(fortfid,"%lf %d %d %lf",&timetemp,&it,&nfil,&defval);
            if (debug>0)printf("Skipping  %d %d %lf %d %d %lf \n",i,i,timetemp,it,nfil,defval); 
	    for (j=0;j<nfil;j++){fscanf(fortfid,"%d %lf %lf",&k,&utemp,&vtemp);}
         }
      }
  
      /* Set outbound arrays */
      if(debug>0)printf("Writing time structure ... \n");	   
      ttime=mxCreateDoubleMatrix(1,ndsets,mxREAL);
      mxFree(mxGetPr(ttime));
      mxSetPr(ttime,time);
      mxSetFieldByNumber(plhs[0], 0, 4, ttime);

      if(debug>0)printf("Writing iter structure ... \n");	   
      iiter=mxCreateDoubleMatrix(1,ndsets,mxREAL);
      mxFree(mxGetPr(iiter));
      mxSetPr(iiter,iter);
      mxSetFieldByNumber(plhs[0], 0, 5, iiter);
      
      if(debug>0)printf("Writing u structure %dx%d ... \n",NN,ndsets);	   
      uu=mxCreateDoubleMatrix(NN,ndsets,mxREAL);
      mxFree(mxGetPr(uu));
      mxSetPr(uu,u);
      mxSetFieldByNumber(plhs[0], 0, 6, uu);
      
      if(debug>0)printf("Writing v structure %dx%d ... \n",NN,ndsets);	   
      vv=mxCreateDoubleMatrix(NN,ndsets,mxREAL);
      mxFree(mxGetPr(vv));
      mxSetPr(vv,v);
      mxSetFieldByNumber(plhs[0], 0, 7, vv);

   } 
   /* --- END CASE FORT.62,64,74 -------------------------------------- */
 
   if (debug>0){
      printf("\nLeaving Mex file...\n");
   }

   fclose(fortfid);
   return;
}
