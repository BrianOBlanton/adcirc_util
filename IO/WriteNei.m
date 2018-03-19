function errno=WriteNei(fgs,nbs)
%WriteNei Write a FEM neighbor file in .nei format.
%   WriteNei requires 2 input arguments
%      1) fgs   - fem_grid_struct
%      2) nbs - node neighbor list
%
%   If fname is omitted, WRITE_NEI enables a file browser
%   with which the user can specify the .NEI file.
%   Otherwise, fname is the name of the .nei file, relative 
%   or absolute (fullpath), including the suffix .'nei'.
%   This input is a string so it must be enclosed in
%   single quotes. 
%
%   WRITE_NEI checks the lengths of input arrays and the
%   structure of the neighbor list for some minimal
%   error checking.  This is not, however, very rigorous.
%
% CALL: err=write_nei(x,y,bc,z,nbs,fname);
%          err=write_nei(x,y,bc,z,nbs);
%
% Written by : Brian O. Blanton 
%              March 1996
%

if nargin==0 & nargout==0
   disp('Call as:  err=WriteNei(fgs,nbs);')
   return
end


if nargin ~=2 
   error(['WriteNei requires 2 input arguments; type "help WRITE_NEI"']);
end

errno=1;

% open fname
[pfid,message]=fopen([fgs.name '.nei'],'w');
if pfid==-1
   error([fpath fname,' not found. ',message]);
end

fprintf(pfid,'%d\n',length(fgs.x));
[nn,nnb]=size(nbs);
fprintf(pfid,'%d\n',nnb);

xmin=min(fgs.x);xmax=max(fgs.x);
ymin=min(fgs.y);ymax=max(fgs.y);
fprintf(pfid,'%.3f   %.3f   %.3f   %.3f\n',[xmax ymax xmin ymin]);

fmt1='%5d %14.6f %14.6f %1d %8.3f ';
fmt2='%5d ';
for i=2:nnb
   fmt2=[fmt2 '%9d '];
end
fmtstr=[fmt1 fmt2 '\n'];
nwrite=nnb+5;
nnn=1:nn;
bc=zeros(size(fgs.x));
fprintf(pfid,eval('fmtstr'),[nnn(:) fgs.x(:) fgs.y(:) bc(:) fgs.z(:) nbs]');

errno=0;
if nargout==0
   clear errno
end

return

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

