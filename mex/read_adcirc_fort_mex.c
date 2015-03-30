/*
 
 read_adcirc_fort_mex.c
 Binary read function for ADCIRC fort files 
 
 the calling function (read_adcirc_fort.m) handles
 most error checking.  
 
 27 Sep, 2005
 
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

/* ---- read_adcirc_fort_mex will be called as :
 D=read_adcirc_fort_mex(fname,unit,debug,stride,strip,level); 
 where fname is the char string containing the fort.?? filename
 and unit is the ?? number.  This selects what gets read.
 ---------------------------- */

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
   const char **fieldnames;
   int *dims;
   
   int i,j,it,k; 
   double *time,*iter;

   char *fortfilename,*header,*linejunk;
   int strlen,status,*supported_units=NULL,nunits,ifortunit;
   double *fortunit=NULL;
   FILE *fortfid=NULL, *fopen();
   double NaN=mxGetNaN();
   double FillThresh=-9999.;
   int SetFill=1;
   int itest;
         
   int NDSETS1,NDSETS2,NSTEP,IFLAG,NN,ndsets;
   int Level=-1,NFEN=1;
   int actual_steps;
   double *DT,*NND,*NFEND;
   mxArray *dt, *nn, *nfen, *ttime, *iiter;
   mxArray *zz,*uu,*vv, *ww;
   double *z,*u,*v, *w, *s;
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
         
   nunits=13;
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
   supported_units[9]=45;
   supported_units[10]=46;
   supported_units[11]=42;
   supported_units[12]=43;
   
   /* --- check I/O arguments ------------------------------------------ */
   if      (nrhs != 8)  mexErrMsgTxt("read_adcirc_fort_mex requires 8 input arguments.");
   else if (nlhs != 1)  mexErrMsgTxt("read_adcirc_fort_mex requires 1 output argument.");
      
   /* --- First input must be a string --------------------------------- */
   if (mxIsChar(prhs[0]) != 1)mexErrMsgTxt("First input must be a string.");

   /* --- Second input must be from supported_units  ------------------- */
   fortunit=mxGetPr(prhs[1]);
   ifortunit=(int)fortunit[0];

   /* --- get debug flag ----------------------------------------------- */
   temp=mxGetPr(prhs[2]);
   debug=(int)temp[0];

   /* --- get stride flag ---------------------------------------------- */
   temp=mxGetPr(prhs[3]);
   Stride=(int)temp[0];

   /* --- get strip flag ----------------------------------------------- */
   temp=mxGetPr(prhs[4]);
   Strip=(int)temp[0];

   /* --- get starting iter number ------------------------------------- */
   temp=mxGetPr(prhs[5]);
   IterStart=(int)temp[0];
   
   /* --- get ending iter number --------------------------------------- */
   temp=mxGetPr(prhs[6]);
   IterEnd=(int)temp[0];
   
   /* get levels to read in, if 3-d file */
   if( (ifortunit == 45) | (ifortunit == 46) ){
      temp=mxGetPr(prhs[7]);
      Level=(int)temp[0]; 
   }

   if (debug>0){
      printf("\nInside Mex file...\n");
      printf("\nnrhs          = %d\n",nrhs);
      printf("Debug         = %d\n",debug);
      printf("Unit number   = %d\n",ifortunit);
      printf("Stride        = %d\n",Stride);
      printf("IterStart     = %d\n",IterStart);
      printf("IterEnd       = %d\n",IterEnd);
      printf("Strip         = %d\n",Strip);
      printf("Level         = %d\n",Level);
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
      return;}
   
   /* --- specify structure field names, ifortunit based --------------- */
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
   else if( (ifortunit == 72) | (ifortunit == 74) ){
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
   else if( (ifortunit == 45) | (ifortunit == 46) ){
      nfields=9;
      fieldnames=mxCalloc(nfields,sizeof(*fieldnames));
      fieldnames[0]="filename";
      fieldnames[1]="header";
      fieldnames[2]="dt";
      fieldnames[3]="nn";
      fieldnames[4]="nfen";
      fieldnames[5]="t";
      fieldnames[6]="uvel";
      fieldnames[7]="vvel";
      fieldnames[8]="wvel";
   }
   else if( (ifortunit == 42) | (ifortunit == 43) ){
      nfields=10;
      fieldnames=mxCalloc(nfields,sizeof(*fieldnames));
      fieldnames[0]="filename";
      fieldnames[1]="header";
      fieldnames[2]="dt";
      fieldnames[3]="nn";
      fieldnames[4]="nfen";
      fieldnames[5]="t";
      fieldnames[6]="uvel";
      fieldnames[7]="vvel";
      fieldnames[8]="wvel";
      fieldnames[9]="z";
   }

   /* --- create return structure -------------------------------------- */
   plhs[0] = mxCreateStructMatrix(1, 1, nfields , fieldnames);
 
   header=mxCalloc(BUFFER_SIZE,sizeof(char));
   linejunk=mxCalloc(BUFFER_SIZE,sizeof(char));
   fgets(header,BUFFER_SIZE,fortfid);  
   if (debug>0)printf("\nfort.XX header line=%s\n",header);

   DT=(double *) mxDvector(0,1);
   NND=(double *) mxDvector(0,1);
   NFEND=(double *) mxDvector(0,1);

   if( ifortunit == 42 | ifortunit == 43 | ifortunit == 45 | ifortunit == 46 ){
                                          /* NDSET3DSV, NSTA3DV, DTDP*NSPO3DSV, NSPO3DSV, NFEN, IRTYPE */
      fscanf(fortfid,"%d %lf %lf %d %lf %d",&NDSETS1,   NND,     DT,            &NSTEP,   NFEND,&IFLAG);}
   else{
      fscanf(fortfid,"%d %lf %lf %d %d",&NDSETS1,NND,DT,&NSTEP,&IFLAG);}

   fgets(linejunk,BUFFER_SIZE,fortfid);  
	   
   /* Set NDSETS2 to min of NDSETS1 and IterEnd */
   NDSETS2=IMIN(NDSETS1,IterEnd);
   
   ndsets=ceil((float)(NDSETS2-Strip)/((float)Stride));

   NN=(int)NND[0];
   NFEN=(int)NFEND[0];

   if (debug>0){
      printf("Total Number of Data Sets in file      = %d\n",  NDSETS1); 
      printf("Max Number of Data Sets to read        = %d\n",  NDSETS2); 
      printf("Data Set Stride                        = %d\n",  Stride);
      printf("Data Set Strip                         = %d\n",  Strip);
      printf("Number of Data Sets to return          = %d\n",  ndsets);
      printf("Number of ADCIRC nodes                 = %d\n",  NN); 
      printf("Number of ADCIRC vertical nodes        = %d\n",  NFEN); 
      printf("Time Step in fort file                 = %lf\n", DT[0]); 
      printf("Number of Model steps per output       = %d\n",  NSTEP);
      printf("Level to Return                        = %d\n",  Level);
      if (SetFill){
         printf("Field with values less than %lf will be set to NaN\n",FillThresh);
      }
   }

   if (Level>NFEN){
      mexErrMsgTxt("Level specified exceeds number of vertical levels in file. Terminal.");
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

   nfen=mxCreateDoubleMatrix(1,1,mxREAL);
   mxFree(mxGetPr(nfen));
   mxSetPr(nfen,NFEND);
	   
   mxSetFieldByNumber(plhs[0], 0, 0, mxCreateString(fortfilename));
   mxSetFieldByNumber(plhs[0], 0, 1, mxCreateString(header));
   mxSetFieldByNumber(plhs[0], 0, 2, dt);
   mxSetFieldByNumber(plhs[0], 0, 3, nn);   

   /* --- switch on fort number ---------------------------------------- */
   
   /* --- Case fort.61,63,69,71,73 ---------------------------------------- */
   if ( ifortunit == 61 | ifortunit == 63 | ifortunit == 69 | ifortunit == 71 | ifortunit == 73 ){      	   
      /* --- allocate return arrays         ---------------------- */
      z=   (double *) mxDvector(0,NN*ndsets);
      time=(double *) mxDvector(0,ndsets);
      iter=(double *) mxDvector(0,ndsets);

      for (j=0;j<ndsets;j++){time[j]=NaN;iter[j]=NaN;}  /* init to NaN */
      for (j=0;j<NN*ndsets;j++){z[j]=NaN;} /* init to NaN */
      
      /* --- Scan data part of fort.?? file ---------------------------- */

      /* Strip off Strip timesteps */
      if (debug>0)printf("Stripping off %d time steps.\n",Strip);
      for (i=0;i<Strip;i++){
         double ztemp;
         fscanf(fortfid,"%lf %d",&timestrip,&it);
	      for (j=0;j<NN;j++){
	         fscanf(fortfid,"%d %lf",&k,&ztemp);
         }
	      if (debug>0)printf("Skipped %d %lf %d %lf\n",i,timestrip,it,ztemp);
      } 

      actual_steps=0;
      for (i=Strip;i<NDSETS2;i++){

         itest=i % Stride;

         if ( itest == 0 ){

            fscanf(fortfid,"%lf %d",&time[actual_steps],&iter[actual_steps]);

            if (debug>1)printf("Scanning %d %d %lf %d\n",
                i,actual_steps,time[actual_steps],iter[actual_steps]); 

            for (j=0;j<NN;j++){fscanf(fortfid,"%d %lf",&k,&z[j+actual_steps*NN]);}

	    if (SetFill){
               for (j=0;j<NN;j++){
                  if (z[j+actual_steps*NN]<FillThresh){
                     z[j+actual_steps*NN]=NaN;
                  } 
	       }
	    }
            actual_steps++;
         }

         /* Skip this iter */
         else{
            double timetemp, ztemp;
            fscanf(fortfid,"%lf %d",&timetemp,&it);
            if (debug>1)printf("Skipped !! %d %lf %d\n",i,timetemp,it); 
            for (j=0;j<NN;j++){fscanf(fortfid,"%d %lf",&k,&ztemp);}
         }
      }

      /* Set outbound arrays */
      if(debug>0)printf("\nWriting time structure ... \n");	   
      ttime=mxCreateDoubleMatrix(1,ndsets,mxREAL);
      mxFree(mxGetPr(ttime));
      mxSetPr(ttime,time);
      mxSetFieldByNumber(plhs[0], 0, 4, ttime);
      
      if(debug>0)printf("\nWriting iter structure ... \n");	   
      iiter=mxCreateDoubleMatrix(1,ndsets,mxREAL);
      mxFree(mxGetPr(iiter));
      mxSetPr(iiter,iter);
      mxSetFieldByNumber(plhs[0], 0, 5, iiter);
      
      if(debug>0)printf("\nWriting zeta structure %dx%d ... \n",NN,ndsets);	   
      zz=mxCreateDoubleMatrix(NN,ndsets,mxREAL);
      mxFree(mxGetPr(zz));
      mxSetPr(zz,z);
      mxSetFieldByNumber(plhs[0], 0, 6, zz);
   
   } 
   /* --- END CASE FORT.61,63,69,71,73 ------------------------------------ */

   /* --- Case fort.62,64,72,74 ---------------------------------------- */
	   
   else if( ifortunit == 62 | ifortunit == 64 | ifortunit == 72 | ifortunit == 74 ){

      /* --- Scan data part of fort.?? file ---------------------------- */
      u=   (double *) mxDvector(0,NN*ndsets);
      v=   (double *) mxDvector(0,NN*ndsets);
      time=(double *) mxDvector(0,ndsets);
      iter=(double *) mxDvector(0,ndsets);

      for (j=0;j<NN*ndsets;j++){u[j]=NaN;} /* init to NaN */
      for (j=0;j<NN*ndsets;j++){v[j]=NaN;} /* init to NaN */
      for (j=0;j<ndsets;j++){time[j]=NaN;iter[j]=NaN;}  /* init to NaN */
	   
      /* Strip off Strip timesteps  */
      for (i=0;i<Strip;i++){
	 double utemp,vtemp;
	 fscanf(fortfid,"%lf %d",&timestrip,&it);
	 for (j=0;j<NN;j++){
	   fscanf(fortfid,"%d %lf %lf",&k,&utemp,&vtemp);
         }
	 if (debug>0)printf("Skipped %d %lf %d %lf %lf\n",i,timestrip,it,utemp,vtemp);
      } 

      actual_steps=0;
      for (i=Strip;i<NDSETS2;i++){
         
         itest=i % Stride;
 
         if ( itest == 0 ){
            fscanf(fortfid,"%lf %lf",&time[actual_steps],&iter[actual_steps]);
            if (debug>1)printf("Reading!! %d %d %lf %d\n",i,actual_steps,time[actual_steps],iter[actual_steps]); 
            for (j=0;j<NN;j++){fscanf(fortfid,"%d %lf %lf",&k,&u[j+actual_steps*NN],&v[j+actual_steps*NN]);}

            for (j=0;j<NN;j++){
               if (u[j+actual_steps*NN]<FillThresh){
		       u[j+actual_steps*NN]=NaN;
		       v[j+actual_steps*NN]=NaN;
	       } 
	    }
            actual_steps++;
         }
	 
 	 /* Skip this iter */
         else{
	    double timetemp, utemp, vtemp;
            fscanf(fortfid,"%lf %d",&timetemp,&it);
	    if (debug>1)printf("Skipped !! %d %lf %d\n",i,timetemp,it); 
	    for (j=0;j<NN;j++){fscanf(fortfid,"%d %lf %lf",&k,&utemp,&vtemp);}
         }
      }
  
      /* Set outbound arrays */
      if(debug>0)printf("\nWriting time structure ... \n");	   
      ttime=mxCreateDoubleMatrix(1,ndsets,mxREAL);
      mxFree(mxGetPr(ttime));
      mxSetPr(ttime,time);
      mxSetFieldByNumber(plhs[0], 0, 4, ttime);

      if(debug>0)printf("Writing iter structure ... \n");	   
      iiter=mxCreateDoubleMatrix(1,ndsets,mxREAL);
      mxFree(mxGetPr(iiter));
      mxSetPr(iiter,iter);
      mxSetFieldByNumber(plhs[0], 0, 5, iiter);

      if(debug>0)printf("\nWriting u structure ... \n");	   
      uu=mxCreateDoubleMatrix(NN,ndsets,mxREAL);
      mxFree(mxGetPr(uu));
      mxSetPr(uu,u);
      mxSetFieldByNumber(plhs[0], 0, 6, uu);

      if(debug>0)printf("\nWriting v structure ... \n");	   
      vv=mxCreateDoubleMatrix(NN,ndsets,mxREAL);
      mxFree(mxGetPr(vv));
      mxSetPr(vv,v);
      mxSetFieldByNumber(plhs[0], 0, 7, vv);

   } 
   /* --- END CASE FORT.62,64,72,74 ------------------------------------ */

   /* --- Case fort.42,43 ---------------------------------------------- */
	   
   else if( ifortunit == 42 | ifortunit == 43 ){
      int dim_s=3*(NFEN-1)+2;
      int dim_uvw=NN*ndsets*NFEN;

      mxSetFieldByNumber(plhs[0], 0, 4, nfen);

      /* --- Scan data part of fort.?? file ---------------------------- */
      s=   (double *) mxDvector(0,dim_s);
      u=   (double *) mxDvector(0,dim_uvw);
      v=   (double *) mxDvector(0,dim_uvw);
      w=   (double *) mxDvector(0,dim_uvw);
      time=(double *) mxDvector(0,ndsets);
      printf("dim_s, dim_uvw= %d %d \n",dim_s,dim_uvw);

      for (j=0;j<dim_uvw;j++){u[j]=NaN;v[j]=NaN;w[j]=NaN;} /* init to NaN */
      for (j=0;j<ndsets;j++){time[j]=NaN;}  /* init to NaN */

      /* Strip off Strip timesteps */
      for (i=0;i<Strip;i++){
	 double utemp, vtemp, wtemp;
         int jj;
	 fscanf(fortfid,"%lf %d",&timestrip,&it);

         /* read in sigma levels */
         for (j=0;j<dim_s;j++){fscanf(fortfid,"%lf",&s[j]);} 

            for (j=0;j<NN;j++){
            fscanf(fortfid,"%d",&k);
            for (jj=0;jj<NFEN;jj++){
                  fscanf(fortfid,"%lf %lf %lf",&utemp, &vtemp, &wtemp);
            }
         }
	 if (debug>1)printf("Skipped : %d %lf %d\n",i,timestrip,it);
      } 

      actual_steps=0;
      for (i=Strip;i<NDSETS2;i++){
         double tempu, tempv, tempw;
         itest=i % Stride;
 
         if ( itest == 0 ){

            if (fscanf(fortfid,"%lf %d",&time[actual_steps],&it)  == EOF ){
		    printf("EOF reached early.\n");
		    break;
	    }
            printf("Reading : %d %lf %d\n",i,time[actual_steps],it); 
            if (debug>1)printf("Reading!! %d %d %lf %d\n",i,actual_steps,time[actual_steps],it); 

            /* read in sigma levels */
            for (j=0;j<dim_s;j++){
               fscanf(fortfid,"%lf",&s[j]); 
            } 
            if (debug>1)printf("%d %lf %lf %lf\n",it,time[actual_steps],s[0],s[dim_s-1]); 

      // int dim_uvw=NN*ndsets*NFEN;
            for (j=0;j<NN;j++){
               fscanf(fortfid,"%d",&k);
               int kk;
               for (kk=0;kk<NFEN;kk++){
           	  int idx;
                  fscanf(fortfid,"%lf %lf %lf",&tempu, &tempv, &tempw);
		  //idx=kk + actual_steps*NN + j*(ndsets*NFEN);
		  idx=actual_steps + j*ndsets + kk*(ndsets*NN);
                  u[idx]=tempu;
                  v[idx]=tempv;
                  w[idx]=tempw;
               }
            }
            actual_steps++;
         }
         else{
	    double timetemp, tempu, tempv, tempw; 
            int jj;

            fscanf(fortfid,"%lf %d",&timetemp,&it);

            /*  read in sigma levels */
            for (j=0;j<dim_s;j++){fscanf(fortfid,"%lf",&s[j]);} 

            for (j=0;j<NN;j++){
               fscanf(fortfid,"%d",&k);
               for (jj=0;jj<NFEN;jj++){
                  fscanf(fortfid,"%lf %lf %lf",&tempu, &tempv, &tempw);
                }
            }
	    if (debug>1)printf("Skipped !! %d %lf %d\n",i,timetemp,it); 
	    printf("Skipped : %d %lf %d\n",i,timetemp,it); 
         }
      }
      
      /* Set outbound arrays */
      if(debug>0)printf("\nWriting time structure ... \n");	   
      ttime=mxCreateDoubleMatrix(1,ndsets,mxREAL);
      mxFree(mxGetPr(ttime));
      mxSetPr(ttime,time);
      mxSetFieldByNumber(plhs[0], 0, 5, ttime);

      dims=mxIvector(0,2);
      dims[0]=NFEN;dims[1]=ndsets;dims[2]=NN; 

      if(debug>0)printf("\nWriting u structure ... \n");	   
      uu=mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      mxFree(mxGetPr(uu));
      mxSetPr(uu,u);
      mxSetFieldByNumber(plhs[0], 0, 6, uu);

      if(debug>0)printf("\nWriting v structure ... \n");	   
      vv=mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      mxFree(mxGetPr(vv));
      mxSetPr(vv,v);
      mxSetFieldByNumber(plhs[0], 0, 7, vv);

      if(debug>0)printf("\nWriting w structure ... \n");	   
      ww=mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      mxFree(mxGetPr(ww));
      mxSetPr(ww,w);
      mxSetFieldByNumber(plhs[0], 0, 8, ww);

      if(debug>0)printf("\nWriting z structure ... \n");	   
      zz=mxCreateDoubleMatrix(1,dim_s,mxREAL);
      mxFree(mxGetPr(zz));
      mxSetPr(zz,s);
      mxSetFieldByNumber(plhs[0], 0, 9, zz);

   }
   /* --- END CASE FORT.42,43 ------------------------------------------ */
   
   /* --- CASE FORT.45,46 ---------------------------------------------- */
	   
   else if( ifortunit == 45 | ifortunit == 46 ){
      int dim_s=3*(NFEN-1)+2;
      int dim_uvw=NN*ndsets*NFEN;

      mxSetFieldByNumber(plhs[0], 0, 4, nfen);

      /* --- Scan data part of fort.?? file ---------------------------- */
      s=   (double *) mxDvector(0,dim_s);
      u=   (double *) mxDvector(0,dim_uvw);
      v=   (double *) mxDvector(0,dim_uvw);
      w=   (double *) mxDvector(0,dim_uvw);
      time=(double *) mxDvector(0,ndsets);
      printf("dim_s, dim_uvw= %d %d \n",dim_s,dim_uvw);

      for (j=0;j<dim_uvw;j++){u[j]=NaN;v[j]=NaN;w[j]=NaN;} /* init to NaN */
      for (j=0;j<ndsets;j++){time[j]=NaN;}  /* init to NaN */

      /* Strip off Strip timesteps */
      for (i=0;i<Strip;i++){
	 double utemp, vtemp, wtemp;
         int jj;
	 fscanf(fortfid,"%lf %d",&timestrip,&it);

         /* read in sigma levels */
         for (j=0;j<dim_s;j++){fscanf(fortfid,"%lf",&s[j]);} 

            for (j=0;j<NN;j++){
            fscanf(fortfid,"%d",&k);
            for (jj=0;jj<NFEN;jj++){
                  fscanf(fortfid,"%lf %lf %lf",&utemp, &vtemp, &wtemp);
            }
         }
	 if (debug>1)printf("Skipped : %d %lf %d\n",i,timestrip,it);
      } 

      actual_steps=0;
      for (i=Strip;i<NDSETS2;i++){
         double *tempu, *tempv, *tempw;
         itest=i % Stride;
 
         if ( itest == 0 ){
            int jj; 

            if (fscanf(fortfid,"%lf %d",&time[actual_steps],&it)  == EOF ){
		    printf("EOF reached early.\n");
		    break;
	    }
            printf("Reading : %d %lf %d\n",i,time[actual_steps],it); 
            if (debug>1)printf("Reading!! %d %d %lf %d\n",i,actual_steps,time[actual_steps],it); 

            /* read in sigma levels */
            for (j=0;j<dim_s;j++){fscanf(fortfid,"%lf",&s[j]);}
            if (debug>1)printf("%d %lf %lf %lf\n",it,time[actual_steps],s[0],s[dim_s-1]); 

            /* double *tempu, *tempv, *tempw; */
            tempu=(double *) mxDvector(0,NFEN);
            tempv=(double *) mxDvector(0,NFEN);
            tempw=(double *) mxDvector(0,NFEN);

            for (j=0;j<NN;j++){
               fscanf(fortfid,"%d",&k);
               for (jj=0;jj<NFEN;jj++){
                  fscanf(fortfid,"%lf %lf %lf",&tempu[jj], &tempv[jj], &tempw[jj]);
               }
               u[j+actual_steps*NN]=tempu[Level-1];
               v[j+actual_steps*NN]=tempv[Level-1];
               w[j+actual_steps*NN]=tempw[Level-1];
            }

            /*for (j=0;j<NN;j++){
               if (u[j+actual_steps*NN]<FillThresh){
		       u[j+actual_steps*NN]=NaN;
		       v[j+actual_steps*NN]=NaN;
		       w[j+actual_steps*NN]=NaN;
	       } 
	    } */

            actual_steps++;
         }
         else{
	    double timetemp;  /*, tempu, tempv, tempw; */
            int jj;

            fscanf(fortfid,"%lf %d",&timetemp,&it);

            /*  read in sigma levels */
            for (j=0;j<dim_s;j++){fscanf(fortfid,"%lf",&s[j]);} 

            for (j=0;j<NN;j++){
               fscanf(fortfid,"%d",&k);
               for (jj=0;jj<NFEN;jj++){
                  fscanf(fortfid,"%lf %lf %lf",&tempu, &tempv, &tempw);
                }
            }
	    if (debug>1)printf("Skipped !! %d %lf %d\n",i,timetemp,it); 
	    printf("Skipped : %d %lf %d\n",i,timetemp,it); 
         }
      }
      
      /* Set outbound arrays */
      if(debug>0)printf("\nWriting time structure ... \n");	   
      ttime=mxCreateDoubleMatrix(1,ndsets,mxREAL);
      mxFree(mxGetPr(ttime));
      mxSetPr(ttime,time);
      mxSetFieldByNumber(plhs[0], 0, 5, ttime);

      if(debug>0)printf("\nWriting u structure ... \n");	   
      uu=mxCreateDoubleMatrix(NN,ndsets,mxREAL);
      mxFree(mxGetPr(uu));
      mxSetPr(uu,u);
      mxSetFieldByNumber(plhs[0], 0, 6, uu);

      if(debug>0)printf("\nWriting v structure ... \n");	   
      vv=mxCreateDoubleMatrix(NN,ndsets,mxREAL);
      mxFree(mxGetPr(vv));
      mxSetPr(vv,v);
      mxSetFieldByNumber(plhs[0], 0, 7, vv);

      if(debug>0)printf("\nWriting w structure ... \n");	   
      ww=mxCreateDoubleMatrix(NN,ndsets,mxREAL);
      mxFree(mxGetPr(ww));
      mxSetPr(ww,w);
      mxSetFieldByNumber(plhs[0], 0, 8, ww);

   }
   /* --- END CASE FORT.45,46 ------------------------------------------ */
   
   if (debug>0){
      printf("\nLeaving Mex file...\n");
   }

   fclose(fortfid);
   return;
}
