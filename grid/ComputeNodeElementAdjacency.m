function [icne, c] = ComputeNodeElementAdjacency(fgs)
%ComputeNodeElementAdjacency compute node/element connectivity
% ComputeNodeElementAdjacency(fem_grid_struct) computes the table of
% elements connected to each node, using the element in the input
% fem_grid_struct. 
%
%  INPUT : fem_grid_struct - (from LOADGRID, ExtractGrid, etc, see FEM_GRID_STRUCT)       
%
% OUTPUT : icne - matrix of size [NN x maxconnections] where maxconnections is
%                 the maximum number of elements that contain any node.
%          c    - number of elements associated with each node ([NN x 1])
%
%   CALL : icne=ComputeNodeElementAdjacency(fgs);
%
% Winter 2022

% check arguments
if nargin ==0 
   disp('icne=ComputeNodeElementAdjacency(fem_grid_struct);')
   return
end  

if ~is_valid_struct(fgs)
   error('Argument to ComputeNodeElementAdjacency must be a valid fem_grid_struct.')
end

icne=zeros(fgs.nn,20);
c=zeros(fgs.nn,1);

for j=1:fgs.ne
%    if rem(i,100)==0, fprintf('%d\n',i); end
    n=fgs.e(j,:);
    % increment counter for nodes n
    c(n)=c(n)+1;
    % insert element j into each node row 
    icne(n(1),c(n(1)))=j;
    icne(n(2),c(n(2)))=j;
    icne(n(3),c(n(3)))=j;
end

icne(:,all(icne==0))=[];
if nargout==1, clear c; end
