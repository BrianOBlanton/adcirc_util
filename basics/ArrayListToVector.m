function out=ArrayListToVector(in)
% ArrayListToVector - convert a  java.util.ArrayList to a MATLAB double
% vector
%
% out=ArrayListToVector(java.util.ArrayList)
%
% Example:
%    A=ArrayList;
%    class(A)
%    for i=1:10
%       add(A,randi(100));
%    end
%    A = [53.0, 78.0, 13.0, 63.0, 35.0, 34.0, 58.0, 87.0, 20.0, 68.0];
%    V=ArrayListToVector(A)'
%    V =
%        53    78    13    63    35    34    58    87    20    68

% 3 different ways of converting.  Number 2 is faster for large lists.

% 1) 

%out=cell2mat(in.toArray.cell);

% 2) 

n = javaArray('java.lang.Double', in.size);
in.toArray(n);
out=double(n);

% 3) 

% out=NaN*ones(in.size,1);
% for i=0:in.size-1
%    out(i+1)=in.get(i);
% end
