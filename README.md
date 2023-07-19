# adcirc_util
## Draft public matlab codes for dealing with ADCIRC

MATLAB functions to contour, plot, etc triangular finite element model output (ADCIRC).

11 Nov 2022

Dependencies: None, but nctoolbox (https://github.com/nctoolbox/nctoolbox) will make life much easier. Nctoolbox greatly simplifies handling ADCIRC netCDF files.

Clone repo and then add <PATH_TO_adcirc_util>/adcirc_util to the MATLAB function path, then execute   

<pre>
  AddAdcircPaths(<PATH_TO_adcirc_util>/adcirc_util)
</pre>

To add ADCIRC paths to the MATLAB MATLABPATH variable in the user's matlab/startup.m file, add:

<pre>
  global ADCIRC
  ADCIRC='<INSTALL_DIR>/adcirc_util';
  AddAdcircPaths(ADCIRC)
</pre>

## Basic commands

### Load a grid
<pre>g=grd_to_opnml('gridname');</pre>

or, using nctoolbox and a url to a ADCIRC netCDF "global" file: 
<pre>
  % locally available maxele file
  url='maxele.63.nc'; 
  
  % or: 
  % maxele file on a TDS
  url='http://tds.renci.org/thredds/dodsC/NCFS_CURRENT_SYNOPTIC/maxele.63.nc';
  
  nc=ncgeodataset(url);
  g=ExtractGrid(nc);
</pre>

### Draw grid elements, add contour
<pre>
  helems=drawelems(g);
  hcont=lcontour(g,'z',0,'LineWidth',1,'color','r');
  axis('equal')
</pre>


