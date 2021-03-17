function fem_grid_struct=grd_to_opnml(fort14name,verbose)
%GRD_TO_OPNML Convert an ADCIRC grd file to an OPNML fem_grid_struct.
% Convert an ADCIRC grd file to an OPNML fem_grid_struct.
% ADCIRC grid information assumed in "fort.14" format.
% The boundary/island information at the tail of the fort.14
% file is ignored.
%
% Input:  fort14name - path/name of fort.14 file;  if not passed,
%                      assumes fort.14 in the currect working dir.
% Output: fem_grid_struct - OPNML grid structure
%
% Call:   fem_grid_struct=grd_to_opnml(fort14name);
%         fem_grid_struct=grd_to_opnml;


if ~exist('verbose','var')
   verbose=false;
end

if ~islogical(verbose)
   error('Verbose arg to grd_to_opnml must be logical')
end

if ~exist('fort14name','var')
   % assume fort.14 filename in the current wd.
   fort14name='fort.14';
end
if verbose, fprintf('Scanning %s ... ',fort14name), end

% Open fort.14 file
[fid,message]=fopen(fort14name,'r');
if (fid<0)
   error(message)
end

% Get grid info
gridname=fgetl(fid);

l=fgetl(fid);
temp=strsplit(strip(l));
ne=str2double(temp{1});
nn=str2double(temp{2});

% Get node locations
if verbose, fprintf('\nnodes = '),end
temp=fscanf(fid,'%d %f %f %f',[4 nn])';
idx=temp(:,1);
x=temp(:,2);
y=temp(:,3);
z=temp(:,4);
if verbose, fprintf('%d ... ',nn),end

% Get elements
if verbose, fprintf('\nelements = '),end 
temp=fscanf(fid,'%d %d %d %d %d',[5 ne])';
e=temp(:,3:5);
if verbose, fprintf('%d ... ',ne),end

if range(e(:)) > length(x)
   etemp=e(:);
    for i=1:length(idx)
        n=idx(i);
        ich=etemp==n;
        etemp(ich)=i;
    end
    e=reshape(etemp,size(e));
end
        
fem_grid_struct.name=strtrim(gridname);
fem_grid_struct.x=x;
fem_grid_struct.y=y;
fem_grid_struct.z=z;
fem_grid_struct.e=e;
fem_grid_struct.bnd=detbndy(e);
fem_grid_struct.nn=length(x);
fem_grid_struct.ne=length(e);

% scan open boundary
if verbose, fprintf('\nNumber of open boundary segments = '), end
fem_grid_struct.nopenboundaries=fscanf(fid,'%d',1);
fgets(fid);

if (fem_grid_struct.nopenboundaries==0)
    fem_grid_struct.nopenboundarynodes{1}=0;
    fem_grid_struct.ob{1}={0};
    fgets(fid);
else    
    fem_grid_struct.elevation=fscanf(fid,'%d',1);fgets(fid);
    for i=1:fem_grid_struct.nopenboundaries
        fem_grid_struct.nopenboundarynodes{i}=fscanf(fid,'%d',1);
        fgets(fid);
        temp=fscanf(fid,'%d',fem_grid_struct.nopenboundarynodes{i});
        fem_grid_struct.ob{i}=temp;
    end
end
if verbose, fprintf('%d ... ',fem_grid_struct.nopenboundaries),end

% scan land boundary
if verbose, fprintf('\nNumber of land boundary segments = '),end 
fem_grid_struct.nland=fscanf(fid,'%d',1);
fgets(fid);
fem_grid_struct.nlandnodes=fscanf(fid,'%d',1);
fgets(fid);

if verbose, fprintf('%d ... ',fem_grid_struct.nland),end
if verbose, fprintf('\n'),end 

n24=0;
n23=0;
%n0=0;
n23nodes=0;
n24pairs=0;

%fem_grid_struct.nfluxnodes=[0];
fem_grid_struct.nlandnodes=0;
fem_grid_struct.ibtype=0;
fem_grid_struct.ln={0};
fem_grid_struct.weirheights={0};
nodeStrings=cell(fem_grid_struct.nland,1);

for i=1:fem_grid_struct.nland
   temp=fscanf(fid,'%d',2);
   rr=fgets(fid); % get remainder of line
   nodeStrings{i}=rr;

   fem_grid_struct.nlandnodes(i)=temp(1);
   fem_grid_struct.ibtype(i)    =temp(2);
   %if verbose, fprintf('\n      %d %d %d ... ',i,temp(1),temp(2)),end

   switch fem_grid_struct.ibtype(i)    % On ibtype
   
   case {0, 1, 2, 10, 11, 12, 20, 21, 22, 30, 52} 
       temp=NaN*ones(fem_grid_struct.nlandnodes(i),1);
       for j=1:fem_grid_struct.nlandnodes(i)
          temp(j)=fscanf(fid,'%d',1);
          fgets(fid);
       end
      %temp=fscanf(fid,'%d',fem_grid_struct.nlandnodes(i));
      fem_grid_struct.ln{i}=temp;

   case {3, 13, 23}          % Exterior Boundary 
      n23=n23+1;
      temp=fscanf(fid,'%d %f %f',[3 fem_grid_struct.nlandnodes(i)])';
      n23nodes=n23nodes+fem_grid_struct.nlandnodes(i);
      fem_grid_struct.ln{i}=temp(:,1);
     
   case {4, 24}          % Node pairs for weirs
      n24=n24+1;
      temp=fscanf(fid,'%d %d %f %f %f',[5 fem_grid_struct.nlandnodes(i)])';
      n24pairs=n24pairs+fem_grid_struct.nlandnodes(i);
      fem_grid_struct.ln{i}=temp(:,1:2);
      fem_grid_struct.weirheights{i}=temp(:,3);
      
   otherwise
      fprintf('Encountered a boundary type (%d) not coded, on segment %d.\n',temp(2),i)
      if temp(2)>100
          fprintf('Boundary type (%d) exceeds allowable range. Terminal.\n',temp(2));
          error('Last Node String: %s',nodeStrings{i-1})
      end
   end
end

if verbose, fprintf('\nNumber of Weir segments = %d \n',n24), end
%if verbose, fprintf('\n Last line of file read: %s',temp);end
fclose(fid);

if isempty(fem_grid_struct.name)
    fem_grid_struct.name='changeme';
end
fem_grid_struct.n23nodes=n23nodes;
fem_grid_struct.n24pairs=n24pairs;
fem_grid_struct.nweir=n24;

try 
    fem_grid_struct=belint(fem_grid_struct);
    fem_grid_struct=el_areas(fem_grid_struct);
catch ME
    fprintf('Returning incomplete fem_grid_struct.\n')
    throw(ME)
end



return
