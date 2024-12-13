function crl=curl(fem_grid_struct,u,v)
%CURL Compute the CURL of a FEM 2-D vector field
% CURL(fem_grid_struct,u,v) computes the curl of a 2-D 
% vector field (u,v) (del X (u,v)) over the FEM domain 
% specified by fem_grid_struct.  The result is a 2-D 
% scalar field containing the curl, returned to the workspace.
%
% The computed quantity is:	crl=dv/dx - du/dy, i.e.,
% the k-th component of (iDx + jDy + kDz) X (iu + jv + kw)
%
%  INPUT : fem_grid_struct (from LOADGRID, see FEM_GRID_STRUCT)
%          u,v - u,v vector field
%
% OUTPUT : crl - del X (u,v)
%
%   CALL : crl=curl(fem_grid_struct,u,v);
%
% ALL ARGUMENTS ARE REQUIRED.
%
% Written by : Brian O. Blanton
% Spring 1998
%

if nargin ~=3 
   error('    CURL MUST (ONLY) HAVE 3 INPUT ARGUMENTS.');
end

if nargout~=1
   error('    CURL must have 1 (and only 1) output argument.')
end

if ~is_valid_struct(fem_grid_struct)
   error('    First argument to CURL must be a valid fem_grid_struct.')
end

crl=curlmex5(fem_grid_struct.x,...
             fem_grid_struct.y,...
	     fem_grid_struct.e,...
	     u,v);
%
%LabSig  Brian O. Blanton
%        Department. of Marine Sciences
%        12-7 Venable Hall
%        CB# 3300
%        University of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        brian_blanton@unc.edu
