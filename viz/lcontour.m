function h=lcontour(fem_grid_struct,Q,cval,varargin)
%LCONTOUR contour a scalar on a FEM grid.
%   LCONTOUR contours a scalar field on the input FEM grid.
%   LCONTOUR accepts a vector of values to be contoured 
%   over the provided mesh.  
%
%   INPUT : fem_grid_struct (from LOADGRID, see FEM_GRID_STRUCT)
%           Q    - scalar to be contoured upon; must be a 1-D vector 
%                  or the single character 'z', IN SINGLE QUOTES!!
%           cval - vector of values to contour
%
%           In order to contour the FEM domain bathymetry, pass
%           in the string 'z' in place of an actual scalar field Q.
%           You could, of course, pass in the actual bathymetry as
%           the scalar to contour.  Otherwise, Q must be a 1-D vector
%           with length equal to the number of nodes in the FEM mesh.
%
%           Any property name/value pair that LINE accepts can be 
%           passed to LCONTOUR. See LINE help for details.
%
%  OUTPUT :  h - the handle to the contour line(s) drawn
%
%    CALL : >> h=lcontour(fem_grid_struct,Q,cval,pn1,pv1,pn2,pv2,...)
%     OR    >> h=lcontour(fem_grid_struct,'z',cval,pn1,pv1,pn2,pv2,...)

% Written by : Brian O. Blanton
% 
% 07 Mar, 2004: moved drawing of contours outside of computational
%               loop to speed up rendering of graphics over slow
%               net connections
% added map catch/linem, 21 Aug 14, BOB

% 
% 
% VERIFY INCOMING STRUCTURE
%
if ~isstruct(fem_grid_struct)
   msg=char(' ',...
               'First argument to LCONTOUR not a structure.  Perhaps its',...
               'the element list.  If so you should use LCONTOUR4, which',...
               'takes the standard grid arrays (e,x,...).  The first ',...
               'argument to LCONTOUR MUST be a fem_grid_struct.',' ');
   disp(msg)
   error(' ')
end
if ~is_valid_struct(fem_grid_struct)
   error('    fem_grid_struct to LCONTOUR invalid.')
end


% Strip off propertyname/value pairs in varargin not related to
% "line" object properties.
k=1;
ax=gca;

while k<length(varargin)
  switch lower(varargin{k})
    case 'axes'
      ax=varargin{k+1};
      varargin([k k+1])=[];
    otherwise
      k=k+2;
  end
end

e=fem_grid_struct.e;
x=fem_grid_struct.x;
y=fem_grid_struct.y;

% DETERMINE SCALAR TO CONTOUR
%
if ischar(Q)
   Q=fem_grid_struct.z;
else
   % columnate Q
   Q=Q(:);
   [nrowQ,~]=size(Q);
   if nrowQ ~= length(x)
      error('Length of scalar must be same length as grid coordinates.');
   end   
end
 
% range of scalar quantity to be contoured; columnate cval
Qmax=max(Q);
Qmin=min(Q);
cval=cval(:);
%h=zeros(size(cval));
ch = true(size(cval));
h = gobjects(size(cval));

for kk=1:length(cval)
%parfor (kk=1:length(cval))
   if (cval(kk) > Qmax) || (cval(kk) < Qmin)
      fprintf('%s not within range of scalar field.  Min = %f  :  Max = %f\n',num2str(cval(kk)),Qmin,Qmax);
        ch(kk)=false;
   else
   
% Call cmex function contmex5
      C=contmex5(x,y,e,Q,cval(kk));
      if(size(C,1)*size(C,2)~=1)
         X = [ C(:,1) C(:,3) NaN*ones(size(C(:,1)))]';
         Y = [ C(:,2) C(:,4) NaN*ones(size(C(:,1)))]';
         XX{kk} = X(:);
         YY{kk} = Y(:);
         %len(kk)=length(X(:));
         %h(kk)=line(XX{kk},YY{kk},'LineStyle','-',varargin{:},'UserData',cval(kk),'Tag','contour','Clipping','on');
         %set(h(kk),'ZData',1*ones(size(get(h(kk),'XData'))))
      else
         disp(['CVal ' num2str(cval(kk)) ' within range but still invalid.']);
         ch(kk)=false;
      end
   end
end 


for kk=1:length(cval)
    if ch(kk)
        if ismap(gca)     
            h(kk)=linem(ax,YY{kk},XX{kk});
        else
            h(kk)=line(ax,XX{kk},YY{kk});
        end
        set(h(kk),'LineStyle','-',varargin{:},'UserData',cval(kk),'Tag','contour','Clipping','on')
        set(h(kk),'DisplayName',sprintf('%.2f',cval))
    end
end
    



% for kk=1:length(cval)
%     if ishandle(h(kk))
%         set(h(kk),'ZData',ones(size(get(h(kk),'XData'))))
%     end
% end
% 
% h(isnan(h))=0;

% 
% if exist('plotfx')
%    plotfx;
% end

return

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
%        Summer 1997
%            Mod 08 Mar, 2004   





 
