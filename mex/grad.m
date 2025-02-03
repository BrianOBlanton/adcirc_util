function gr=grad(e,x,y,q)
%GRAD Compute the gradient of a FEM 2-D scalar field
% GRAD(e,x,y,s) computes the gradient of a 2-D 
% scalar field (s) over the FEM domain specified by the 
% element list (e) and the corresponding horizontal node
% coordinates (x,y).  The result is a 2-D vector field
% returned to the workspace. 
%
% Call as: gd=grad(e,x,y,s);
%
% All arguments are REQUIRED.
%

if nargout~=1
   error('GRAD must have 1 (and only 1) output argument.')
end

if (size(x)~=size(y))|(size(x)~=size(q))
   error('Size of x,y, and q MUST be equal.');
end

[grdu grdv]=gradmex5(x,y,e,q);
gr=grdu(:)+sqrt(-1)*grdv(:);


%
%LabSig  Brian O. Blanton
%        Department of Marine Sciences
%        12-7 Venable Hall
%        CB# 3300
%        University of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        brian_blanton@unc.edu
%
