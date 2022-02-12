function [dst,inearest]=nearest_node(fgs,x,y,n)
%NEAREST_NODE find FEM node nearest input coordinates
%
% Input:   fgs - fem_grid_struct
%            x - east coordinates
%            y - north coordinates
%            n - number of neighbors
% Output: dst  - distance to nearest node
%         inearest - nearest node
% Call as: [dst,inearest]=nearest_node(fgs,x,y,n);

dst=NaN*ones(length(x),n);
inearest=dst;

for i=1:length(x)
   xdif=fgs.x-x(i);
   ydif=fgs.y-y(i);
   [a,b]=sort(xdif.*xdif+ydif.*ydif);
%    idx=a==0;
%    a(idx)=[];
%    b(idx)=[];
   dst(i,:)=a(1:n);
   inearest(i,:)=b(1:n);
end


