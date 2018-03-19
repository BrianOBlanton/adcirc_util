function retval=colormesh3d(fem_grid_struct,Q,sarg)

% COLORMESH3D draw FEM mesh in 3-d given, colored by a scalar field.
%  
%     INPUT : fem_grid_struct (from LOADGRID, see FEM_GRID_STRUCT)
%             Q         - scalar to color with (optional)
%             "nofill"  - flag (optional)
%
%             With no scalar specified to contour, COLORMESH2D
%             defaults to the bathymetry fem_grid_struct.z
%
%   Output : hp - vector of handles, one for each element patch drawn.
%
%            COLORMESH3D colors the mesh using the scalar Q.  If Q
%            is omitted from the argument list, COLORMESH3D uses
%            the bathymetry z to color the surface.  
%
%            The default viewpoint is Azimuth = -27. degrees, 
%            Elevation = 30 degrees.  See the MATLAB VIEW command for 
%            more information.
%
%  Call as : >> hp=colormesh3d(fgs,Q,'nofill')
%
%   Author : Brian O. Blanton
%

% DEFINE ERROR STRINGS
err1=['more than 1 column in x-coordinate vector'];
err2=['more than 1 column in y-coordinate vector'];
err3=['more than 1 column in z-coordinate vector'];
err4=['node coordinate vectors must be the same length (len(x)~=len(y))'];
err5=['node coordinate vectors must be the same length (len(x|y)~=len(z))'];
err6=['scalar to be contoured must be 1-D'];
err7=['length of scalar must be the same length as coordinate vectors'];
err8=str2mat('??? Error using ==> colormesh3d',...
             'COLORMESH3D is confused by the number of columns in ',...
             'the element file.  It must be only three (i1 i2 i3) ',...
             'or four (node# i1 i2 i3)');
err9=str2mat('??? Error using ==> colormesh3d',...
             'COLORMESH3D needs atleast 4 input arguments.',...
             'Input : elems - element list (.ele or .tri type)',...
             '        x     - x-coordinate list',...
             '        y     - y-coordinate list',...
             '        z     - z-coordinate list (bathymetry)',...
             '        Q     - scalar to color with (optional)',...
             '    "nofill"  - no-interior triangles flag (optional)');
err10=str2mat('??? Error using ==> colormesh3d',...
             'COLORMESH3D optional flag argument must',...
             'be the string "nofill".');
             
             
% VERIFY INCOMING STRUCTURE
%
if ~isstruct(fem_grid_struct)
   error(err11)
end
if ~is_valid_struct(fem_grid_struct)
   error('    fem_grid_struct to COLORMESH2D invalid.')
end
 
e=fem_grid_struct.e;
x=fem_grid_struct.x;
y=fem_grid_struct.y;
z=fem_grid_struct.z;

x=x(:);
y=y(:);
z=z(:);

if exist('Q')
   if isstr(Q)
      if strcmp(Q,'nofill')
         Q=z;
         sarg='nofill';
      else
         disp(err10)
         return
      end
   else
      Q=Q(:);
   end
else
   Q=z;
end

[nrowQ,ncolQ]=size(Q);
if ncolQ~=1 && nrowQ~=1
   error(err6);
end
if nrowQ ~= length(x)
   error(err7)
end

e=e';

% eliminate elements with neg areas
iding=find(fem_grid_struct.ar<0);
if ~isempty(iding)
    fprintf('%s','There are negative areas in the grid.  These are being eliminated ',...
        'from the patch matrix, assuming they are elements that span a branch cut.')
    e(:,iding)=[];
end

[m,n]=size(e);
xt=x(e);
yt=y(e);
zt=z(e);
Qt=Q(e);
xt=reshape(xt,m,n);
yt=reshape(yt,m,n);
zt=reshape(zt,m,n);
Qt=reshape(Qt,m,n);

% delete previous colorsurf objects
% The commented line below will only work in 4.2c or greater.
%delete(findobj(gca,'Type','patch','Tag','colorsurf'))

  
% Create patch object
%
if exist('sarg','var')
   if strcmp(sarg,'nofill')
      hp=patch(xt,yt,zt,Qt,'EdgeColor','interp',...
               'FaceColor','none','Tag','colorsurf');
   end
else
   hp=patch(xt,yt,zt,Qt,'EdgeColor','interp','Tag','colorsurf');
end


% Output if requested.
if nargout==1,retval=hp;end
%
%        Brian O. Blanton
%        Curriculum in Marine Science
%        Ocean Processes Numerical Modeling Laboratory
%        15-1A Venable Hall
%        CB# 3300
%        Uni. of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        919-962-4466
%        blanton@marine.unc.edu
%
%        Jan 1995
%
