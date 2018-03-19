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


f = netcdf.create(f_filename, 'WRITE');
f.title                   = f_title;
f.institution             = f_institution;
f.institution_code        = f_institution_code;
f.institution_OPeNDAP_url = f_institution_OPeNDAP_url;
f.contact                 = f_contact;
f.project                 = f_project;
f.project_url             = f_project_url;
f.model_domain            = fem_grid_struct.name;
f.Conventions             = 'CF-x.x' ;
f.format_category         = 'finite element triangular grid' ;
f.grid_type               = 'Triangular' ;


% dimensions

f('node') = length(fem_grid_struct.x);
f('nele') = size(fem_grid_struct.e,1);
f('nface')=size(fem_grid_struct.e,2);
f('nbd')=length(fem_grid_struct.bnd);
f('nbi')=4;
f('charlen')=132;

% variables
f{'lon'}=ncdouble('node');
  f{'lon'}.long_name='Longitude';
  f{'lon'}.units='degrees_east';
  f{'lon'}.standard_name='longitude';
  f{'lon'}.axis='X';
f{'lat'}=ncdouble('node');
  f{'lat'}.long_name='Latitude';
  f{'lat'}.units='degrees_north';
  f{'lat'}.standard_name='latitude';
  f{'lat'}.axis='Y';
f{'depth'}=ncdouble('node');
  f{'depth'}.long_name='Bathymetry';
  f{'depth'}.units='meters';
  f{'depth'}.standard_name='depth';
  f{'depth'}.grid='tri_grid';
%   f{'depth'}.axis='Z';

f{'tri_grid'}=ncchar('charlen');
  f{'tri_grid'}.domain_name=fem_grid_struct.name;
  f{'tri_grid'}.grid_name='triangular_mesh';
  f{'tri_grid'}.Horizontal_Triangular_Element_Incidence_List='ele';
  f{'tri_grid'}.Boundary_Segment_Node_List = 'bnd' ;
  f{'tri_grid'}.Index_start='1';
  f{'tri_grid'}.grid_type='Triangular';

f{'ele'}=ncint('nele','nface');
  f{'ele'}.long_name='Horizontal_Triangular_Element_Incidence_List';
  f{'ele'}.standard_name='XXX';

f{'bnd'}=ncint('nbd','nbi');
  f{'bnd'}.long_name='Boundary_Segment_Node_List';
  f{'bnd'}.standard_name='XXX';

f{'land_binary_mask'}=ncdouble('node');
  f{'land_binary_mask'}.long_name='land_binary_mask';
  f{'land_binary_mask'}.standard_name='land_binary_mask';

f{'water_binary_mask'}=ncdouble('node');
  f{'water_binary_mask'}.long_name='water_binary_mask';
  f{'water_binary_mask'}.standard_name='water_binary_mask';

f{'lon'}(:)    =fem_grid_struct.x;
f{'lat'}(:)    =fem_grid_struct.y;
f{'depth'}(:)  =fem_grid_struct.z;
f{'ele'}(:)    =fem_grid_struct.e;
f{'bnd'}(:,1:2)=fem_grid_struct.bnd;

fillvals=ncfillvalues;
delete ncfillvalues.nc

idx=1*(fem_grid_struct.z>0);
idx(idx==0)=fillvals.ncdouble;
f{'water_binary_mask'}(:)=idx;   
idx=1*(fem_grid_struct.z<=0);
idx(idx==0)=fillvals.ncdouble;
f{'land_binary_mask'}(:)=idx;   

% close object
f = close(f);

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
%        WINTER 2005
%

