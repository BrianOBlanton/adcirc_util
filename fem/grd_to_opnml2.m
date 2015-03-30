function fem_grid_struct=grd_to_opnml2(fort14name,verbose)
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
% Call:   fem_grid_struct=grd_to_opnml2(fort14name);
%         fem_grid_struct=grd_to_opnml2;


if ~exist('verbose')
   verbose=true;
end

if ~islogical(verbose)
   error('Verbose arg to grd_to_opnml must be logical')
end

verboseI=0;
if verbose,verboseI=1;,end


if ~exist('fort14name')
   % assume fort.14 filename in the current wd.
   fort14name='fort.14';
end
if verbose, fprintf('Scanning %s : \n',fort14name), end

% test fort.14 file
if ~exist(fort14name,'file')
   error('Can not find %s file',fort14name)
end

fem_grid_struct=read_adcirc_grd_mex(fort14name,verboseI);

fem_grid_struct.name=strtrim(fem_grid_struct.name);
if verbose, fprintf('Adding boundary segments ...\n'),end 
fem_grid_struct.bnd=detbndy(fem_grid_struct.e);
fem_grid_struct=belint(fem_grid_struct);
if verbose, fprintf('Adding element areas ...\n'),end 
fem_grid_struct=el_areas(fem_grid_struct);
if verbose, fprintf('Adding element centroids ...\n'),end 
fem_grid_struct=attach_elem_centroids(fem_grid_struct);
%if verbose, fprintf('Adding STRtree ...\n'),end 
%fem_grid_struct.strtree=ComputeStrTree(fem_grid_struct);
%if verbose, fprintf('Adding Element Adjacency ...\n'),end 
%fem_grid_struct.icee=ComputeElementAdjacency(fem_grid_struct);
return

% scan open boundary
if verbose, fprintf('   open boundary = '), end
fem_grid_struct.nopen=fscanf(f14,'%d',1);fgets(f14);
fem_grid_struct.elevation=fscanf(f14,'%d',1);fgets(f14);
if (fem_grid_struct.nopen==0)   
   fem_grid_struct.nopennodes={0};
   fem_grid_struct.ob={0};  
else

for i=1:fem_grid_struct.nopen
       fem_grid_struct.nopennodes{i}=fscanf(f14,'%d',1);
       fgets(f14);
       temp=fscanf(f14,'%d',fem_grid_struct.nopennodes{i});
       fem_grid_struct.ob{i}=temp;
    end
end
if verbose, fprintf('%d ... ',fem_grid_struct.nopen),end

% scan land boundary
if verbose, fprintf('\nland boundary segments = '),end 
fem_grid_struct.nland=fscanf(f14,'%d',1);fgets(f14);
fem_grid_struct.nlandnodestotal=fscanf(f14,'%d',1);fgets(f14);
if verbose, fprintf('%d ... ',fem_grid_struct.nland),end

n24=0;
n23=0;
n0=0;
n23nodes=0;
n24pairs=0;

fem_grid_struct.nlandnodes=[0];
fem_grid_struct.ibtype=[0];
fem_grid_struct.ln={0};
fem_grid_struct.weirheights={0};

for i=1:fem_grid_struct.nland

   temp=fscanf(f14,'%d',2);
   fgets(f14); % get remainder of line

   fem_grid_struct.nlandnodes(i)=temp(1);
   fem_grid_struct.ibtype(i)    =temp(2);

   switch temp(2)   % On ibtype
   
   case {0, 1, 2, 10, 11, 12, 20, 21, 22, 30, 52} 
      temp=fscanf(f14,'%d',fem_grid_struct.nlandnodes(i));
      fem_grid_struct.ln{i}=temp;
      
   case {3, 13, 23}          % Exterior Boundary 
      n23=n23+1;
      temp=fscanf(f14,'%d %f %f',[3 fem_grid_struct.nlandnodes(i)])';
      n23nodes=n23nodes+fem_grid_struct.nlandnodes(i);
      fem_grid_struct.ln{i}=temp(:,1);
     
   case {4, 24}          % Node pairs for weirs
      n24=n24+1;
      temp=fscanf(f14,'%d %d %f %f %f',[5 fem_grid_struct.nlandnodes(i)])';
      n24pairs=n24pairs+fem_grid_struct.nlandnodes(i);
      fem_grid_struct.ln{i}=temp(:,1:2);
      fem_grid_struct.weirheights{i}=temp(:,3);
      
   otherwise
      disp(['Boundary type not coded: ' int2str(temp(2))])
   end
end

fclose(f14);

fem_grid_struct.n23nodes=n23nodes;
fem_grid_struct.n24pairs=n24pairs;
fem_grid_struct.nweir=n24;

fem_grid_struct=belint(fem_grid_struct);
fem_grid_struct=el_areas(fem_grid_struct);


if verbose, fprintf('Number of Weir segments = %d \n',n24), end

return

