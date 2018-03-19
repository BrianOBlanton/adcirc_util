function rv1=colormesh2de(fem_grid_struct,Q,varargin)
%COLORMESH2DE draw an element-based FEM mesh in 2-d colored by a scalar quantity.
%
%   INPUT : fem_grid_struct (from LOADGRID, see FEM_GRID_STRUCT)
%	    Q	      - scalar to color with, [1 X ne] (required)
%
%  OUTPUT : hp - vector of handles to patchs drawn.
%
%  PN/PV pairs accepted by COLORMESH2D
%	    NBand     - number of contour bands to compute  (def=16);
%           ColorMap  - colormap to use (def=jet(nband);)
%  All other PVPN pairs are passed to PATCH.
%
%   CALL : >> hp=colormesh2de(fem_grid_struct,Q,pn1,pv1,...);
%
% Written by : Brian O. Blanton
%

narginchk(1,2);

% VERIFY INCOMING STRUCTURE
%
if ~isstruct(fem_grid_struct)
   error('First argument to COLORMESH2DE must be a structure.')
end
if ~is_valid_struct(fem_grid_struct)
   error('fem_grid_struct to COLORMESH2DE invalid.')
end
 
e=fem_grid_struct.e;
x=fem_grid_struct.x;
y=fem_grid_struct.y;

NBand=16;
ColorMap='jet';
FaceColor='flat';
EdgeColor='none';

[nrowQ,~]=size(Q);
if nrowQ ~= length(e)
   error('length of scalar must be the same length as element list')
end
Q=Q(:);

% Strip off parameter/value pairs in varargin not related to
% "line" object properties.
k=1;
while k<length(varargin),
   switch lower(varargin{k}),
   case 'nband',
      NBand=varargin{k+1};
      varargin([k k+1])=[];
   case 'colormap',
      ColorMap=varargin{k+1};
      varargin([k k+1])=[];
   case 'facecolor',
      FaceColor=varargin{k+1};
      varargin([k k+1])=[];
   case 'edgecolor',
      EdgeColor=varargin{k+1};
      varargin([k k+1])=[];
   otherwise
     k=k+2;
   end;
end;
if length(varargin)<2
   varargin={};
end

% delete previous colorsurf objects
%delete(findobj(gca,'Type','patch','Tag','colorsurf'))

hp=patch('faces',e,'vertices',[x y],'facevertexcdata',Q,'Tag','colorsurf','FaceColor',FaceColor,...
'EdgeColor',EdgeColor,varargin{:});
eval(['colormap(' ColorMap  '(' sprintf('%d',NBand)  ')'  ')' ])


% Output if requested.
if nargout==1,rv1=hp;end

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

 
