function metis_diag_error(ErrorFile,RunDir)

if ~exist('RunDir','var')
   RunDir='./';
end
%disp(RunDir)
if ~exist(RunDir,'dir')
   error('RunDir DNE. Terminal.')
end

cd(RunDir)

fid=fopen(ErrorFile);
fclose(fid);





