function [Conn,kk]=bnd_conn(fem_struct)
%
% BND_CONN computes ordered boundary connectivity lists 
%	from the fem_struct.bnd as computed by DETBNDY.
%	The structure returned has KK fields named conn.bndi,
%	where i=1:kk.  By default, the boundary with largest area
%	is named the land boundary, listed CCW, and given first.
%	Subsequent boundaries are island boundaries and listed CW.
%
% 	BND_CONN requires the following inputs:
%		fem_struct
%
%	output vars:
%		conn	- structure containing ordered boundary
%			  lists for each closed boundary
%		kk	- number of boundaries found
% 
%	[conn,kk]=bnd_conn(fem_struct);
%
% Calls: coastline/comparea
%
% Catherine R. Edwards
% Last modified: 31 Jul 2001
%
	
bnd=fem_struct.bnd;
tmpbnd=bnd;
kk=0;
%maxar=0;
while(~isempty(tmpbnd))
  connected=0;
  strnn=tmpbnd(1,1); 
  nextnn=tmpbnd(1,2);
  %conbnd=NaN*ones(length(fem_struct.bnd),1);
  conbnd=[tmpbnd(1,1) nextnn]; 
  tmpbnd=tmpbnd(2:end,:);
  %prevnn=0;

    while ~connected

       [i,j]=find(tmpbnd==nextnn);
       if length(i)>1
           jj=rem(j,2)+1;
           for ii=1:length(i)
              possible_connections(ii)=tmpbnd(i(ii),jj(ii));
           end
           dx=fem_struct.x(possible_connections)-fem_struct.x(nextnn);
           dy=fem_struct.y(possible_connections)-fem_struct.y(nextnn);
           ang=atan2(dy,dx)*180/pi;
           idx=ang<0;
           ang(idx)=ang(idx)+360;
           [minang,idx]=min(ang);
           i=i(idx);
       end

       connr=setdiff(tmpbnd(i,:),nextnn);
       %disp([prevnn nextnn connr])
       conbnd=[conbnd connr];
       tmpbnd(i,:)=[];
       if connr==strnn
          connected=1; 
          kk=kk+1;
          %line(fem_struct.x(conbnd),fem_struct.y(conbnd),'LineStyle','-','Color',col(kk),'LineWidth',2)
          %drawnow
       end
       prevnn=nextnn;
       nextnn=connr;
    end
    
    % do area check to determine whether points are ordered CW/CCW; flip CW
    conbnd=conbnd(:);
    ar=comparea(fem_struct.x(conbnd),fem_struct.y(conbnd)); 
    iscw=sign(ar);
    if(iscw>0)
        conbnd=flipud(conbnd);
    end
  
    Conn{kk}=conbnd;

end	



function area = comparea(x,y,dim)
%COMPAREA Area of polygon.
%
%   COMPAREA is a modified version of the Matlab supplied POLYAREA. The
%      area of a polygon is negative if the points are ordered
%      counter-clockwise

if nargin==1, error('Not enough inputs.'); end

if ~isequal(size(x),size(y)), error('X and Y must be the same size.'); end

if nargin==2,
  [x,nshifts] = shiftdim(x);
  y = shiftdim(y);
elseif nargin==3,
  perm = [dim:max(length(size(x)),dim) 1:dim-1];
  x = permute(x,perm);
  y = permute(y,perm);
end

siz = size(x);
if ~isempty(x),
  area = reshape((sum( (x([2:siz(1) 1],:) - x(:,:)).* ...
                 (y([2:siz(1) 1],:) + y(:,:)))/2),[1 siz(2:end)]);
else
  area = sum(x); % SUM produces the right value for all empty cases
end

if nargin==2,
  area = shiftdim(area,-nshifts);
elseif nargin==3,
  area = ipermute(area,perm);
end
