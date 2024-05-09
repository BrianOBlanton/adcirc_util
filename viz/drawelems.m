function hel=drawelems(fem_grid_struct,varargin)
%DRAWELEMS draw 2-D FEM element configuration 
%
% DRAWELEMS draws element boundries given a valid grid structure.  
%
%  INPUT : fem_grid_struct - (from LOADGRID, see FEM_GRID_STRUCT)       
%           
% OUTPUT : hel - handle to the element object.
%
%   CALL : hel=drawelems(fem_grid_struct,p1,v1,...);
%
% Summer 1997
%     

%REARTH=6367500;
%LabelElems='no';

% DEFINE ERROR STRINGS
%err1=['Not enough input arguments; need a fem_grid_struct'];

% check arguments
if nargin ==0 
   disp('h=drawelems(fem_grid_struct,p1,v1,...);')
   return
end  

if ~is_valid_struct(fem_grid_struct)
   error('    Argument to DRAWELEMS must be a valid fem_grid_struct.')
end

% Default propertyname values
TheseElems=[];
MeshHeight=1.;
SPH=true;
ax=[];

% Strip off propertyname/value pairs in varargin not related to
% "line" object properties.
k=1;
while k<length(varargin)
  switch lower(varargin{k})
    case 'meshheight'
      MeshHeight=varargin{k+1};
      varargin([k k+1])=[];
    case 'sph'
      SPH=varargin{k+1};
      varargin([k k+1])=[];
%     case 'labelelems'
%       LabelElems=varargin{k+1};
%       varargin([k k+1])=[];
    case 'theseelems'
      TheseElems=varargin{k+1};
      varargin([k k+1])=[];
    case 'axes'
      ax=varargin{k+1};
      varargin([k k+1])=[];
    otherwise
      k=k+2;
  end
end
if isempty(ax)
    ax=gca;
end

if ~SPH && ~isfield(fem_grid_struct,'xecen_cart')
    error('SPH set to false, but cartesian coord fields are not in input fem_grid_struct.')
end

if length(varargin)<2
   varargin={};
end

% Extract grid fields from fem_grid_struct
elems=fem_grid_struct.e;
if SPH
    x=fem_grid_struct.x;
    y=fem_grid_struct.y;
else
    x=fem_grid_struct.x_cart;
    y=fem_grid_struct.y_cart;
end
z=fem_grid_struct.z;

% eliminate elements outside of current view; if element centroids are
% attached to grid struct
if isempty(TheseElems)
    if isfield(fem_grid_struct,'xecen')
        axx=axis(ax);
        dx=axx(2)-axx(1);
        dy=axx(4)-axx(3);
        ikeep=fem_grid_struct.xecen>axx(1)-dx & ...
              fem_grid_struct.xecen<axx(2)+dx & ...
              fem_grid_struct.yecen>axx(3)-dy & ...
              fem_grid_struct.yecen<axx(4)+dy;
        if ~all(ikeep==0)
            elems=elems(ikeep,:);
        end
    end
else
    elems=elems(TheseElems,:);
end

% COPY FIRST COLUMN TO LAST TO CLOSE ELEMENTS
edges=[elems(:,[1 2]); elems(:,[2 3]) ; elems(:,[3 1])];
[m,~]=size(edges);
xt=[x(edges) NaN*ones([m 1])];
yt=[y(edges) NaN*ones([m 1])];
zt=ones(size([z(edges) NaN*ones([m 1])]));

xt=xt';
yt=yt';
zt=zt';

xt=xt(:);
yt=yt(:);
zt=zt(:);

if ismap(gca)
    mstruct=gcm;
    [xt,yt] = projfwd(mstruct,yt,xt);
end

% DRAW GRID
%hel=line(ax,xt,yt,zt,'Color','k',varargin{:},'Tag','elements');
hel=line(ax,xt,yt,'Color','k',varargin{:},'Tag','elements');



% Default fontsize=20;
% ps=20;
% if strcmp(LabelElems,'yes') && ~isempty(TheseElems)
%    % Build string matrix
%    strlist=num2str(TheseElems,10);
%    % label only those nodes that lie within viewing window.
%    htext=text(x,y,strlist,...
%                     'FontSize',ps,...
%                     'HorizontalAlignment','center',...
%                     'VerticalAlignment','middle',...
%                     'Color','k',...
%                     varargin{:},...
%                     'Tag','Node #');
% end




%hel=line(xt,yt,zt,varargin{:},'Tag','elements');
%set(hel,'ApplicationData',zt);

% if (abs(MeshHeight)>eps)
%    nz=MeshHeight*ones(size(get(hel,'XData')));
%    set(hel,'ZData',nz)
% end


% if exist('plotfx')
%    plotfx;
% end




% if SPH
%     
%     % eliminate edges that cross the branch cut
% %     temp=xt(:,1).*xt(:,2);
% %     iding=find(temp<-1000);
% %     xt(iding,:)=[];
% %     yt(iding,:)=[];
%     
%     
%     xt=xt';
%     yt=yt';
%     zt=zt';
%     
%     xt=xt(:);
%     yt=yt(:);
%     zt=zt(:);
% else
%     xt=xt';
%     yt=yt';
%     %zt=zt';
%     
%     xt=xt(:);
%     yt=yt(:);
%     %zt=zt(:);
%     %[x,y,z]=sph2cart(REARTH*xt*pi/180,REARTH*yt*pi/180,REARTH);
%     [x,y,z]=sph2cart(xt*pi/180,yt*pi/180,1);
%     hel=line(x,y,z,'Color','k',varargin{:},'Tag','elements');
% end







%
%Sig  Brian O. Blanton
%     Renaissance Computing Institute
%     University of North Carolina
%     Chapel Hill, NC
%     brian_blanton@renci.org
%
