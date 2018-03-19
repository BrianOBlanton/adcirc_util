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

/* PROTOTYPES */
/* void BDS_unpack(float *flt, unsigned char *bits, unsigned char *bitmap,
        int n_bits, int n, double ref, double scale);
*/

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
   int *dims;
   
   double *DT;
   mxArray *dt;
   
   int nfields;
   nfields=1;
   fieldnames[0]="dt";

   dims=mxIvector(0,1);
   dims[0]=1;
   dt=mxCreateNumericArray(1,dims,mxDOUBLE_CLASS,mxREAL);
   mxSetPr(dt,DT);

   plhs[0] = mxCreateStructMatrix(1, 1, nfields , fieldnames);
   mxSetFieldByNumber(plhs[0], 0, 0, dt);
   
   return;

}
  
