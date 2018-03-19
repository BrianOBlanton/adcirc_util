function retval=plotbndadc(fem_grid_struct,varargin)
%PLOTBNDADC plot color-coded boundary of FEM mesh
%   PLOTBNDADC(fem_grid_struct) draws a color-coded boundary of 
%   an ADCIRC FEM domain.  
%
%  
%   INPUT : fem_grid_struct (from LOADGRID, see FEM_GRID_STRUCT)
%           standard LINE parameter/value pairs can be passed in.
%
%  OUTPUT : hboun - handle to boundary object drawn
%
%    CALL : hboun=plotbnd(fem_grid_struct,p1,v1,...)
%
% Written by : Brian O. Blanton
% Summer 2007
%

% added varargin, 28 Jun 02, BOB

%if nargin ~=1 
%   error('    Incorrect number of input arguments to PLOTBND');
%end

if ~is_valid_struct(fem_grid_struct)
   error('    Argument to PLOTBNDADC must be a valid fem_grid_struct.')
end

% Extract grid fields from fem_grid_struct
%
bnd=fem_grid_struct.bnd;
x=fem_grid_struct.x;
y=fem_grid_struct.y;

% this is the complete boundary
ns=bnd(:,1);
ne=bnd(:,2);
X=[x(ns) x(ne) NaN*ones(size(ns))]';
Y=[y(ns) y(ne) NaN*ones(size(ns))]';
X=X(:);
Y=Y(:);
hboun=line(X,Y,'Tag','boundary','Color','k',varargin{:});

% legstr={'Land'};
% hleg=hboun;

legstr={};
hleg=[];

% color bnd types
if isfield(fem_grid_struct,'nopenboundaries')
   if fem_grid_struct.nopenboundaries>0
      for i=1:fem_grid_struct.nopenboundaries
         iob=fem_grid_struct.ob{i};
         hob(i)=line(fem_grid_struct.x(iob), fem_grid_struct.y(iob),'Color','b','LineStyle','-','Marker','*');
      end
      legstr={'Elevation'};
      hleg=[hleg;hob(1)];
   end
end

cols={
      'r' 'none' '.'   % type 0
      'm' 'none' '.'   % type 1      
      '-' 'none' '-'   % type 2
      '-' 'none' '-'   % type 3
      '-' 'none' '-'   % type 4
      '-' 'none' '-'   % type 5
      '-' 'none' '-'   % type 6 (n/a)
      '-' 'none' '-'   % type 7 (n/a)
      '-' 'none' '-'   % type 8 (n/a)
      '-' 'none' '-'   % type 9 (n/a)
      'b' 'none' '.'   % type 10 
      'c' 'none' '.'   % type 11
      'c' 'none' '.'   % type 12
      'c' 'none' '.'   % type 13
      '-' 'none' '-'   % type 14 (n/a)
      '-' 'none' '-'   % type 15 (n/a)
      '-' 'none' '-'   % type 16 (n/a)
      '-' 'none' '-'   % type 17 (n/a)
      '-' 'none' '-'   % type 18 (n/a)
      '-' 'none' '-'   % type 19 (n/a)
      'c' 'none' '.'   % type 20
      'c' 'none' '.'   % type 21
      'c' 'none' '.'   % type 22
      'c' 'none' '.'   % type 23
      };

legs={'Land (no normal, free tangential)'
      'Island (no normal, free tangential)'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      'Land (no normal, no tangential)'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      'River (weak normal, free tangential)'
    };

if isfield(fem_grid_struct,'ibtype')

% Boundary types (not open) in grid
Ibtypes=unique(fem_grid_struct.ibtype)

for i=1:length(Ibtypes)
   disp(sprintf(' Drawing ibtype=%d',Ibtypes(i)))
   
   idx=find(fem_grid_struct.ibtype==Ibtypes(i));
   size(idx)
   switch Ibtypes(i)

      case 0  % land (no normal, free tangential)

         j=fem_grid_struct.ln{idx(i)};
         x=fem_grid_struct.x(j);
         y=fem_grid_struct.y(j);
         h=line(x,y,'Color',cols{Ibtypes(i)+1,1},'LineStyle',cols{Ibtypes(i)+1,2},'Marker',cols{Ibtypes(i)+1,3});
         hleg=[hleg;h];
         legstr={legstr{:},legs{i}};
         
      case 1  % island

         for j=1:length(idx)
            k=fem_grid_struct.ln{idx(j)};
            x=fem_grid_struct.x(k);
            y=fem_grid_struct.y(k);
            h=line(x,y,'Color',cols{Ibtypes(i)+1,1},'LineStyle',cols{Ibtypes(i)+1,2},'Marker',cols{Ibtypes(i)+1,3});
         end
         hleg=[hleg;h];
         legstr={legstr{:},'Island'};
              
       case 10  % land (no normal, no tangential)

         for j=1:length(idx)
            k=fem_grid_struct.ln{idx(j)};
            x=fem_grid_struct.x(k);
            y=fem_grid_struct.y(k);
            h=line(x,y,'Color',cols{Ibtypes(i)+1,1},'LineStyle',cols{Ibtypes(i)+1,2},'Marker',cols{Ibtypes(i)+1,3});
         end
         hleg=[hleg;h];
         legstr={legstr{:},'Land (no normal, no tangential)'};
         
      case 22  % "river"

         for j=1:length(idx)
            k=fem_grid_struct.ln{idx(j)};
            x=fem_grid_struct.x(k);
            y=fem_grid_struct.y(k);
            h=line(x,y,'Color',cols{Ibtypes(i)+1,1},'LineStyle',cols{Ibtypes(i)+1,2},'Marker',cols{Ibtypes(i)+1,3});
            h=line(x,y,'Color','c','LineStyle','none','Marker','*');
         end
         hleg=[hleg;h];
         legstr={legstr{:},'River'};
         
      case 24
         i24=find(fem_grid_struct.ibtype==24);
         for ii=1:length(i24)
            j=fem_grid_struct.ln{i24(ii)};
            x=fem_grid_struct.x(j);
            y=fem_grid_struct.y(j);
            line(x,y,'Color','r')

         %   j=floor(length(x)/2);
         %   text(x(j),y(j),int2str(i),'Horiz','center','Vert','middle')

            j=fem_grid_struct.ln{i24(i)};
            x=fem_grid_struct.x(j);
            y=fem_grid_struct.y(j);
            line(x,y,'Color','b')
         end

      otherwise
   end
   
   
end
end

if ~isempty(legstr)
   legend(hleg,legstr)
end

if nargout==1,retval=hboun;,end
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
%        Summer 1997
%
