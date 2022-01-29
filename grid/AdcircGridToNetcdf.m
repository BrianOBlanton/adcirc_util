function AdcircGridToNetcdf(fem_grid_struct,varargin)
%AdcircGridToNetcdf package fem grid in netCDF, with "CF" conventions.
% AdcircGridToNetcdf packages a fem_grid_struct into a 
% netCDF formatted file with loose-CF conVentions.
%
% The following global attributes can be set by passing in 
% the parameter with a relevant value:
% filename      resulting netCDF filename 
%               def=[fem_grid_struct.name '.nc'];
% title         tital global attribute
%               def=['Finite element COMPUTATIONAL GRID netCDF
%               FILE: ' fem_grid_struct.name];
% institution             def='not specified';
% institution_code        def='not specified';
% institution_OPeNDAP_url def='not specified';
% contact                 def='not specified';
% project                 def='not specified';
% project_url             def='not specified';
%
% Call as: AdcircGridToNetcdf(fem_grid_struct,varargin)

if ~is_valid_struct(fem_grid_struct)
    error('    Argument to AdcircGridToNetcdf must be a valid fem_grid_struct.')
end

%default values
[pth,name,ext]=fileparts(fem_grid_struct.name);
f_filename=[name '.nc'];
f_title=['Finite element COMPUTATIONAL GRID netCDF FILE: ' fem_grid_struct.name];
f_institution             ='not specified';
f_institution_code        ='not specified';
f_institution_OPeNDAP_url ='not specified' ;
f_contact                 ='not specified';
f_project                 ='not specified';
f_project_url             ='not specified';

% parse varargin
k=1;
while k<length(varargin)
  switch lower(varargin{k})
     case 'filename'
        f_filename=varargin{k+1};
        varargin([k k+1])=[];
     case 'title'
        f_title=varargin{k+1};
        varargin([k k+1])=[];
     case 'institution'
        f_institution=varargin{k+1};
        varargin([k k+1])=[];
     case 'institution_code'
        f_institution_code=varargin{k+1};
        varargin([k k+1])=[];
     case 'institution_opendap_url'
        f_institution_OPeNDAP_url=varargin{k+1};
        varargin([k k+1])=[];
     case 'project'
        f_project=varargin{k+1};
        varargin([k k+1])=[];
     case 'project_url'
        f_project_url=varargin{k+1};
        varargin([k k+1])=[];
     case 'contact'
        f_contact=varargin{k+1};
        varargin([k k+1])=[];
  otherwise
     k=k+2;   
  end
end

if ~ strcmp(f_filename(end-1:end),'nc') 
    f_filename=[f_filename '.nc'];
end

f = netcdf.create(f_filename, 'WRITE');
GlobalId=netcdf.getConstant('GLOBAL');
netcdf.putAtt(f,GlobalId,'creation_date',datestr(now));
netcdf.putAtt(f,GlobalId,'title',f_title);
netcdf.putAtt(f,GlobalId,'model_domain',strrep(fem_grid_struct.name,'.grd',''));

% f.title                   = f_title;
% f.institution             = f_institution;
% f.institution_code        = f_institution_code;
% f.institution_OPeNDAP_url = f_institution_OPeNDAP_url;
% f.contact                 = f_contact;
% f.project                 = f_project;
% f.project_url             = f_project_url;
% f.model_domain            = fem_grid_struct.name;
% f.Conventions             = 'CF-x.x' ;
% f.format_category         = 'finite element triangular grid' ;
% f.grid_type               = 'Triangular' ;


% dimensions
nn=size(fem_grid_struct.x,1);
ne=size(fem_grid_struct.e,1);
nf=size(fem_grid_struct.e,2);
nb=size(fem_grid_struct.bnd,1);
nbi=4;
DimNod=netcdf.defDim(f,'node',nn);
DimEle=netcdf.defDim(f,'nele',ne);
DimFac=netcdf.defDim(f,'nvertex',nf);
DimNbd=netcdf.defDim(f,'nbd',nb);
DimNbi=netcdf.defDim(f,'nbi',nbi);

% variables
VarLon = netcdf.defVar(f,'x','double',DimNod);
         netcdf.putAtt(f,VarLon, 'long_name',     'longitude');
         netcdf.putAtt(f,VarLon, 'units',         'degrees_east');
         netcdf.putAtt(f,VarLon, 'standard_name', 'longitude');
         netcdf.putAtt(f,VarLon, 'axis',          'X');
         netcdf.putAtt(f,VarLon, 'location',      'node');
         netcdf.putAtt(f,VarLon, 'positive',      'east');

VarLat = netcdf.defVar(f,'y','double',DimNod);
         netcdf.putAtt(f,VarLat, 'long_name',     'latitude');
         netcdf.putAtt(f,VarLat, 'units',         'degrees_north');
         netcdf.putAtt(f,VarLat, 'standard_name', 'latitude');
         netcdf.putAtt(f,VarLat, 'axis',          'Y'); 
         netcdf.putAtt(f,VarLat, 'location',      'node');
         netcdf.putAtt(f,VarLat, 'positive',      'north');

VarDep = netcdf.defVar(f,'depth','double',DimNod);
         netcdf.putAtt(f,VarDep, 'long_name',     'distance below geoid');
         netcdf.putAtt(f,VarDep, 'units',         'm');
         netcdf.putAtt(f,VarDep, 'standard_name', 'depth below geoid');
         netcdf.putAtt(f,VarDep, 'axis',          'Z'); 
         netcdf.putAtt(f,VarDep, 'coordinates',   'time y x');
         netcdf.putAtt(f,VarDep, 'location',      'node');

    %netcdf.putAtt(f,VarDep,'grid','tri_grid'); 

% f{'tri_grid'}=ncchar('charlen');
%   f{'tri_grid'}.domain_name=fem_grid_struct.name;
%   f{'tri_grid'}.grid_name='triangular_mesh';
%   f{'tri_grid'}.Horizontal_Triangular_Element_Incidence_List='ele';
%   f{'tri_grid'}.Boundary_Segment_Node_List = 'bnd' ;
%   f{'tri_grid'}.Index_start='1';
%   f{'tri_grid'}.grid_type='Triangular';

VarEle = netcdf.defVar(f,'ele','int',[DimEle DimFac]);
netcdf.putAtt(f,VarEle,'long_name','Horizontal_Triangular_Element_Incidence_List');

% f{'ele'}=ncint('nele','nface');
%   f{'ele'}.long_name='Horizontal_Triangular_Element_Incidence_List';
%   f{'ele'}.standard_name='XXX';
% 
VarBnd = netcdf.defVar(f,'bnd','double',[DimNbd DimNbi]);
netcdf.putAtt(f,VarBnd,'long_name','Boundary_Segment_Node_List');

% f{'bnd'}=ncint('nbd','nbi');
%   f{'bnd'}.long_name='Boundary_Segment_Node_List';
%   f{'bnd'}.standard_name='XXX';

% f{'land_binary_mask'}=ncdouble('node');
%   f{'land_binary_mask'}.long_name='land_binary_mask';
%   f{'land_binary_mask'}.standard_name='land_binary_mask';
% 
% f{'water_binary_mask'}=ncdouble('node');
%   f{'water_binary_mask'}.long_name='water_binary_mask';
%   f{'water_binary_mask'}.standard_name='water_binary_mask';

netcdf.endDef(f);

netcdf.putVar(f,VarLon,0,nn,fem_grid_struct.x);
netcdf.putVar(f,VarLat,0,nn,fem_grid_struct.y);
netcdf.putVar(f,VarDep,0,nn,fem_grid_struct.z);
netcdf.putVar(f,VarEle,[0 0],[ne nf],fem_grid_struct.e);
netcdf.putVar(f,VarBnd,[0 0],[nb 2],fem_grid_struct.bnd);
 

% fillvals=ncfillvalues;
% delete ncfillvalues.nc

% idx=1*(fem_grid_struct.z>0);
% idx(idx==0)=fillvals.ncdouble;
% f{'water_binary_mask'}(:)=idx;   
% idx=1*(fem_grid_struct.z<=0);
% idx(idx==0)=fillvals.ncdouble;
% f{'land_binary_mask'}(:)=idx;   

% close object
netcdf.close(f);

return

%
%        Brian O. Blanton
%        RENCI
%        University of North Carolina
%        Chapel Hill, NC
%
%        brian_blanton@renci.org
%
%        March 2018

