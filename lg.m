function g=lg(dpath)
%LG autoload a grid in the current working directory
% LG loads a fem_grid_struct in specified directory.  If more than
% one fem_grid_struct exists in the directory, a menu is generated to
% select the grid name from.
%
% Inputs:  directory_path - optional, def='./'
%
% Outputs: fem_grid_struct - ...
%
% Call as: g=lg(dpath);

% Written by: Brian Blanton, Spring '03

%%

if ~exist('dpath')
   dpath='./';
else
   if ~(exist(dpath)==7)
      error('Input directory path to LG is not a directory.')
   end
   dpath=[dpath '/'];
end

dpath

g=[];

% get docking state
dks=get(0,'DefaultFigureWindowStyle'); 
set(0,'DefaultFigureWindowStyle','normal');

nod=dir([dpath '*.nod']);
grd=dir([dpath '*.grd']);
f14=dir([dpath '*.14']);
if isempty(nod) && isempty(grd) && isempty(f14)
   disp(['No grids in ' dpath])
   [griddirs,grids,archdir]=listgrids;
   filelist=strcat(griddirs,strcat('/', grids));
   filelist{end+1}='CANCEL';
   K = menu(['Choose a grid from ' archdir],filelist);
   if strcmp(filelist{K},'CANCEL')
        return
   end
   name=[archdir '/' filelist{K}];

else
   nod=[nod; grd; f14];
   K=1;
   if length(nod)<2
      [~,name,ext] =fileparts(nod(1).name);
      if strcmp(ext,'.grd')
	 name=[name '.grd'];
      elseif strcmp(ext,'.14')
	 name=[name '.14'];
      end
   else
%        for i=1:length(nod)
%            [pathstr,filelist{i},ext] =fileparts(nod(i).name);
%            if strcmp(ext,'.grd')
%                filelist{i}=[filelist{i} '.grd'];
%            end
%        end
%        filelist{end+1}='CANCEL';
       filelist={nod.name 'CANCEL'};
       K = menu('Choose a grid:',filelist);
       if strcmp(filelist{K},'CANCEL')
           return
       end
       name=[dpath filelist{K}];
   end
end

[~,~,ext] =fileparts(name);
%disp(['Loading ' name])
if strcmp(ext,'.grd') || strcmp(ext,'.14') 
   g=grd_to_opnml(name,true);
else
   g=loadgrid(name);
end


g=attach_elem_centroids(g);

set(0,'DefaultFigureWindowStyle',dks);
