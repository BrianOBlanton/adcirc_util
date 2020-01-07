function rv1=colormesh2d(fem_grid_struct,Q,nband)
%COLORMESH2D draw a FEM mesh in 2-d colored by a scalar quantity.
%
%   INPUT : fem_grid_struct (from LOADGRID, see FEM_GRID_STRUCT)
%	    Q	      - scalar to color with (optional)
%	    nband     - number of contour bands to compute (optional)
%
%	    With no scalar specified to contour, COLORMESH2D
%	    defaults to the bathymetry fem_grid_struct.z
%
%  OUTPUT : hp - vector of handles, one for each element patch drawn.
%
%   COLORMESH2D colors the mesh using the scalar Q.  If Q
%   is omitted from the argument list, COLORMESH2D draws
%   the element connectivity in black and white.
%
%   CALL : >> hp=colormesh2d(fem_grid_struct,Q,nband)
%
% Written by : Brian O. Blanton
% added map catch/linem, 21 Aug 14, BOB
%


narginchk(1,3);

% VERIFY INCOMING STRUCTURE
%
if ~isstruct(fem_grid_struct)
   error('First argument to COLORMESH2D must be a structure.')
end
if ~is_valid_struct(fem_grid_struct)
   error('fem_grid_struct to COLORMESH2D invalid.')
end
 
e=fem_grid_struct.e;
x=fem_grid_struct.x;
y=fem_grid_struct.y;

% DETERMINE SCALAR TO CONTOUR
%
if ~exist('Q','var')
  Q=fem_grid_struct.z;
  nband=16;
elseif ischar(Q)
  if strcmpi(Q,'z')
    Q=fem_grid_struct.z;            % Default to bathymetry
  else
     error('Second arg to COLORMESH2D must be ''z'' for depth')
  end
  nband=16;
elseif length(Q)==1
  % nband pass in as Q
  nband=Q;
  Q=fem_grid_struct.z;
else
   % columnate Q
   Q=Q(:);
   [nrowQ,~]=size(Q);
   if nrowQ ~= length(x)
      error('Length of scalar must equal number of nodes in grid.');
   end 
   if nargin==2,nband=16;end
end

if nargin==3
   if length(nband)>1
      error('nband argument to COLORMESH2D must be 1 integer')
   end
end
            
[nrowQ,ncolQ]=size(Q);
if ncolQ>1,error(err4);end
if nrowQ ~= length(x)
   error('length of scalar must be the same length as coordinate vectors')

end
Q=Q(:);

% delete previous colorsurf objects
delete(findobj(gca,'Type','patch','Tag','colorsurf'))

z=0*ones(size(x));

try 
    mm=gcm; %#ok<NASGU>
%    disp('map')
    %hp=patchesm('faces',e,'vertices',[y x z],'facevertexcdata',Q,'EdgeColor','none',...
    %         'FaceColor','interp','Tag','colorsurf');
    hp=patchm('XData',y,'YData',x,'ZData',z,'CData',Q,'Tag','colorsurf');

catch
%   disp('plain')
    hp=patch('XData',x,'YData',y,'ZData',z,'faces',e,'CData',Q, 'EdgeColor','none',...
             'FaceColor','interp','Tag','colorsurf');
     
end

%dcm_obj = datacursormode(gcf);
%set(dcm_obj,'UpdateFcn',{@myupdatefcn,Q,fem_grid_struct})

%colormap(jet(nband))

% Output if requested.
if nargout==1,rv1=hp;end

%if exist('plotfx')
%    plotfx;
%end


function txt = myupdatefcn(~,event_obj,Q,g)
% Customizes text of data tips
pos = get(event_obj,'Position');
%I = get(event_obj, 'DataIndex');
[~,I]= min(abs((g.x-pos(1)).^2+(g.y-pos(2)).^2));
txt = {['lon: ',num2str(pos(1))],...
       ['lat: ',num2str(pos(2))],...
       ['node: ',int2str(I)],...
       ['val: ',num2str(Q(I))]};

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

 
