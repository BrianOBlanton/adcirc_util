function D=read_adcirc_fort20(fname,nb,starttime)
% D=read_adcirc_fort20(fname,nb);
% nb = number of nodes in flux boundary
% starttime = datenum time of start of file, can be empty

if ~exist('starttime','var')
    starttime=0;
elseif isempty(starttime)
    starttime=0;
end

fid=fopen(fname,'r');
if fid < 0
    error('Count not open %s file',fname)
end

D.fname=fname;
D.nb=nb;

D.dt=fscanf(fid,'%f',1);
temp=fscanf(fid,'%f');
nt= length(temp)/nb;
D.z=reshape(temp,[nb nt])';

D.t=(0:D.dt:D.dt*(nt-1))/86400+starttime;
