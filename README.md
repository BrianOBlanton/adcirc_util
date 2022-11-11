# adcirc_util
Draft public matlab codes for dealing with ADCIRC 

Github Repository of MATLAB utility routines for the model ADCIRC.
11 Nov 2022

matlab   - MATLAB functions to contour, plot, etc triangular finite element model output (ADCIRC).
           Add <INSTALL_DIR>/adcirc_util/matlab to the MATLAB function path  

%%% MATLAB %%%
Compile the mex files in adcirc_util/mex.  From within MATLAB, cd to the mex directory and type 
>> makemex 

To add ADCIRC paths to the MATLAB MATLABPATH variable in the user's matlab/startup.m file, add:

% Add distribition paths

global ADCIRC
ADCIRC='<INSTALL_DIR>/adcirc_util';
AddAdcircPaths(ADCIRC)



