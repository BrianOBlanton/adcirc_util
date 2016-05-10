function phi=basis2d(fgs,xylist,j)
%BASIS2D compute basis functions for input points in FEM grid
%   BASIS2D computes the FEM basis functions for a given 
%   horizontal position, specified either in the argument
%   list or with the mouse.
%
%   In determining which element has been selected, 
%   BASIS2D needs elemental areas and shape functions.
%   Element areas are returned by LOADGRID.
%   The routine BELINT computes shape function information
%   and attaches it to the input fem_grid_struct.
%   These two functions MUST be run before BASIS2D will 
%   run.
%
%   BELINT is run as: 
%      new_struct=belint(fem_grid_struct); 
%   If ever needed, EL_AREAS is run as:
%      [new_struct,ineg]=el_areas(fem_grid_struct); 
%
%   INPUT : fem_grid_struct - (from LOADGRID, see FEM_GRID_STRUCT)     
%           xylist  (op)    - points to find elements for [n x 2 double] 
%           j       (op)    - element list corresponding to the xylist
%                             set of points.  Optional, but if passed in, 
%                             points will not be relocated. length(j) must
%                             equal length(xylist).
%           
%   OUTPUT : basis function(s) and element number(s)
%            If elements were NOT passed in, then providing two output
%            arguments will collect the elements determined to contain the 
%            specified points.
% 
%   CALL : >> [phi,j]=basis2d(fem_grid_struct)   for interactive
%    or
%          >> [phi,j]=basis2d(fem_grid_struct,xylist)        
%    or
%          >> phi=basis2d(fem_grid_struct,xylist,j)        
%

% Written by : Brian O. Blanton 
% Summer 1998


if nargin==0 & nargout==0
   disp('phi=basis2d2(fem_grid_struct,j)')
   disp('OR: ')
   disp('[phi,j]=basis2d2(fem_grid_struct)')
   return
end

% Input arguemnt number check
nargchk(1,3,nargin);

if nargin==3 & nargout==2
   error('cannot input AND output element list to BASIS2D')
end

if nargin==1
elseif nargin==2 | nargin==3
end

ifin=find(~isnan(j));

if isempty(ifin)
   phi=[];
   j=[];
   return;
end

xp=xylist(:,1);
yp=xylist(:,2);

phi=NaN*ones(length(j),3);

% Extract local information
n3=fgs.e(j(ifin),:);
x=fgs.x(n3);
if length(xp)==1,x=x';,end
x1=x(:,1);x2=x(:,2);x3=x(:,3);
y=fgs.y(n3);
if length(xp)==1,y=y';,end
y1=y(:,1);y2=y(:,2);y3=y(:,3);
area=fgs.ar(j(ifin));

xptemp=xp(ifin);
yptemp=yp(ifin);

% Basis function #1
a=(x2.*y3-x3.*y2)./(2.0*area);
b=(y2-y3)./(2.0*area);
c=-(x2-x3)./(2.0*area);
phi(ifin,1)=a+b.*xptemp+c.*yptemp;

% Basis function #2
a=(x3.*y1-x1.*y3)./(2.0*area);
b=(y3-y1)./(2.0*area);
c=-(x3-x1)./(2.0*area);
phi(ifin,2)=a+b.*xptemp+c.*yptemp;

% Basis function #3
a=(x1.*y2-x2.*y1)./(2.0*area);
b=(y1-y2)./(2.0*area);
c=-(x1-x2)./(2.0*area);
phi(ifin,3)=a+b.*xptemp+c.*yptemp;

if nargout==0
   clear j phi
end

%
%  Brian O. Blanton
%  Renaissance Computing Institute
%  University of North Carolina
%  Chapel Hill, NC
%
%  Brian_Blanton@Renci.Org
%
%  orig: Summer 1998

