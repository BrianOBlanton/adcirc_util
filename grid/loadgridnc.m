function fgs=loadgridnc(gridname)
%LOADGRIDNC 
%
%  Input: gridname - name of domain to load grid files.  Assumes .nc file
%         extension if none is specified.
%
% Output: LOADGRIDNC returns a fem_grid_struct containing (atleast)
%         the following fields:
%         1)  .name   name of FEM domain
%         2)  .e      node connectivity array (cols 2-4 of .ele
%         file )
%         3)  .x      x-horizontal node list  (col 2 of .nod file )
%         4)  .y      y-horizontal node list  (col 3 of .nod file )
%         5)  .z      bathymetry list         (col 2 of .bat) file
%         6)  .bnd    boundary node/pair list (.bnd file )
%         7)  .nn     number of horizontal nodes
%         8)  .ne     number of horizontal elements
%
%   #7,8 are only for completeness.



global DOMAIN GRIDS GRIDDIRS
         
if nargout ~=1
   error('One output argument required.');
end
         

if nargin==0

   gridname=uigetfile({'*.nc;*.cdf;*.netcdf'},'Click on nc filename');
   if gridname==0
      fgs=[];
      return
   end
  
end



[fpath,name,ext,versn]=fileparts(gridname);
if isempty(ext),ext='.nc';,end

% if absolute or relative path, fpath nonempty, load and exit
if exist(gridname)
   disp(sprintf('%s exists.  Loading ...',gridname))
   fgs=loadgnc(gridname);
   return
else
   % check in griddirs
   imatch=strmatch(name,GRIDS);
   if imatch==0
   
      error(sprintf('No match for %s found in grid database.',name))
     
   end
   
   fname=[deblank(GRIDDIRS(imatch,:)) '/' name   ext];
   fgs=loadgnc(fname);
end



% private functions
function fgs=loadgnc(f_filename)


f = netcdf(f_filename, 'nowrite');
fgs.name=f.model_domain(:);
fgs.x=f{'lon'}(:);
fgs.y=f{'lat'}(:);
fgs.z=f{'depth'}(:);
fgs.e=f{'ele'}(:);
fgs.bnd=f{'bnd'}(:);
fgs.nn=length(fgs.x);
fgs.ne=length(fgs.e);


return

