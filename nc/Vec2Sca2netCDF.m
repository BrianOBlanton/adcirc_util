function Vec2Sca2netCDF(fem_grid_struct,nc,varargin) 
%Var2netCDF - writes a (time-dependent) array to a netCDF "fort.63.nc" file
% 
%  INPUT : fem_grid_struct (from LOADGRID, ExtractGrid, etc. See FEM_GRID_STRUCT)  
%          nc -ncgeodataset of vector field
%          
% OUTPUT : file to disk
%
% PN/PV pairs accepted by Var2netCDF:
%    varname  - variable name in netCDF file to write (def='need_a_better_name')
%    units    - units for variable (def = 'need_better_units';
%    filename - netCDF file name (def=[varname '.nc']).  User must include
%               the ".63.nc" in "filename"
%
%   CALL : Var2netCDF(fem_grid_struct,nc,p1,v1,...)
%

% Winter 2022
% Written by : Brian O. Blanton

if ~is_valid_struct(fem_grid_struct)
    error('Argument to Var2netCDF must be a valid fem_grid_struct.')
end

% parse varargin
varname='need_a_better_name';
units='need_better_units';
filename=[varname '.nc'];
k=1;
while k<length(varargin)
    switch lower(varargin{k})
        case 'varname'
            varname=varargin{k+1};
            varargin([k k+1])=[];
        case 'filename'
            filename=varargin{k+1};
            varargin([k k+1])=[];
        case 'units'
            units=varargin{k+1};
            varargin([k k+1])=[];
        otherwise
            k=k+2;
    end
end

if ~isa(nc,'ncgeodataset')
    error('Input nc variable must be of class ncgeodataset.')
end

% extract time from ncgeodataset
t=nctime(nc);

AdcircGridToNetcdf(fem_grid_struct,'filename',filename);

% open netcdf file in append mode
ncid = netcdf.open(filename,'WRITE');
% set to definition mode
netcdf.reDef(ncid)

% define time as unlimited dim
DimTim=netcdf.defDim(ncid,'time',netcdf.getConstant('UNLIMITED'));

DimNod=netcdf.inqDimID(ncid,'node');
[~,NN]=netcdf.inqDim(ncid,DimNod);

VarID = netcdf.defVar(ncid,varname, 'double', [DimNod DimTim]);
    netcdf.putAtt(ncid,VarID, 'long_name',     varname);
    netcdf.putAtt(ncid,VarID, 'units',         units);
    netcdf.putAtt(ncid,VarID, 'standard_name', varname);
    netcdf.putAtt(ncid,VarID, 'location',      'node');
    netcdf.putAtt(ncid,VarID, 'coordinates',   'time y x');

%    netcdf.putAtt(ncid,VarID,'axis','Y'); 
%    netcdf.putAtt(ncid,VarID,'grid','tri_grid'); 

TimID = netcdf.defVar(ncid, 'time', 'double', DimTim);
    netcdf.putAtt(ncid,TimID, 'long_name', 'model time');
    unitsstr=['seconds since ' datestr(t(1))]; 
    netcdf.putAtt(ncid,TimID, 'units', unitsstr);
    netcdf.putAtt(ncid,TimID, 'standard_name', 'time');
    netcdf.putAtt(ncid,TimID, 'base_date', datestr(t(1),'YYYY-MM-DD HH:mm:ss'));
    %netcdf.putAtt(ncid,TimID, 'axis', 'T'); 

netcdf.endDef(ncid)
u=nc{'windx'};
v=nc{'windy'};

% put data along unlimited dim 
tsecs=seconds(t-t(1));
for i=1:length(tsecs)
    disp(['writing ' int2str(i)] )
    netcdf.putVar(ncid,TimID,i-1,1,tsecs(i));
    
    s=abs(u(i,:)+1i*v(i,:));
    
    netcdf.putVar(ncid,VarID,[0 i-1],[NN 1],s);
    
end

netcdf.close(ncid);
