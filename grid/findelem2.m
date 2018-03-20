function jret=findelem2(fem_grid_struct,xylist,tolerance)
%FINDELEM2 find elements from element hash table
%   FINDELEM2 finds the element number arbitrary locations,
%   using a lookup table to narrow the element search.
%   NaN is returned for each point outside of the FEM domain.
%
%   To determine elements FINDELEM2 needs elemental areas,
%   shape functions, and a hash table for element lookups.
%   The routines BELINT and EL_AREAS compute these arrays
%   and add them to a previously created fem_grid_struct.
%   LOADGRID attaches these fields upon loading a grid. 
%   It also attaches an existing hash table for a given grid,
%   
%   To attach directly:
%     new_struct=belint(fem_grid_struct);
%     [new_struct,ineg]=el_areas(fem_grid_struct);
%     new_struct=gen_element_hash(fem_grid_struct,10);
%
%   INPUT : fem_grid_struct - (from LOADGRID, see FEM_GRID_STRUCT)
%           xylist          - points to find elements for [n x 2 double]
%
%   OUTPUT : an element number(s)
% 
%   CALL : >> j=findelem2(fem_grid_struct)   for interactive
%     OR   >> j=findelem2(fem_grid_struct,xylist)        
%     OR   >> j=findelem2(fem_grid_struct,xylist,tolerance)        
%

%   Written by : Brian O. Blanton 
%   Fall 2005
%

if nargin==0
   disp('Call as: j=findelem2(fem_grid_struct,xylist,[tolerance]);')
   return
end

% VERIFY INCOMING STRUCTURE
%
if ~is_valid_struct(fem_grid_struct)
   error('    fem_grid_struct to FINDELEM invalid.')
end

% Make sure additional needed fields of the fem_grid_struct
% have been filled.
if ~is_valid_struct2(fem_grid_struct)
   error('    fem_grid_struct to FINDELEM invalid.')
end

if ~isfield(fem_grid_struct,'ehash')
   disp('No element hash table attached. Use FINDELEM or GEN_ELEMENT_HASH.')
   return
end
ehash=fem_grid_struct.ehash;

% default tolerance for basis function evaluation
DefaultTolerance=1.e-6;
if ~exist('tolerance')
   tolerance=DefaultTolerance;
end
   

if exist('xylist')
   xp=xylist(:,1);
   yp=xylist(:,2);
%   line(xp,yp,'LineStyle','.','Marker','+')
else
   disp('Click on element ...');
   waitforbuttonpress;
   Pt=gcp;
   xp=Pt(2);yp=Pt(4);
   line(xp,yp,'LineStyle','+')
end

%indices into hash table
jp=ceil((xp-ehash.x(1))/ehash.dx2);
ip=ceil((yp-ehash.y(1))/ehash.dy2);

ifind=find(ip>0 & jp>0);
xtemp=xp(ifind);
ytemp=yp(ifind);
jret=NaN*ones(size(xp));
itemp=ip(ifind);
jtemp=jp(ifind);

ii=unique(itemp);
jj=unique(jtemp);

for j=1:length(jj)
   for i=1:length(ii)
      iii=ii(i);
      jjj=jj(j);
      idx=find(itemp==iii & jtemp==jjj);
      xx=xtemp(idx);
      yy=ytemp(idx);
      %line(xx,yy,'Marker','.','Color','r','LineStyle','none','MarkerSize',1)
      %drawnow
      jsearch=ehash.e{iii,jjj};
      %disp(sprintf('%d %d %d',iii,jjj,length(jsearch)))
      j4=findelemex52(...
               xx,yy,...
               fem_grid_struct.ar,...
               fem_grid_struct.A,...
               fem_grid_struct.B,...
               fem_grid_struct.T,...
               jsearch,tolerance);
      jret(ifind(idx))=j4;
   end
end

return

