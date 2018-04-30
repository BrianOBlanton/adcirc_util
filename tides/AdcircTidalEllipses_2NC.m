function AdcircTidalEllipses_2NC(fgs,D)

% append the grid like this: ncks -A ec2012_v3d.nc f54.nc

% output a netcdf file
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLOBBER'));
f = netcdf.create('tidal_ellipse_parameters.nc', mode);
GlobalId=netcdf.getConstant('GLOBAL');

netcdf.putAtt(f,GlobalId,'creation_date',datestr(now));
netcdf.putAtt(f,GlobalId,'title','ADCIRC ec2015 East Coast Tidal Database');
netcdf.putAtt(f,GlobalId,'model_domain',strrep(fgs.name,'.grd',''));
netcdf.putAtt(f,GlobalId,'Conventions','CF-1.6, UGRID-1.0');
netcdf.putAtt(f,GlobalId,'cdm_data_type','ugrid');
netcdf.putAtt(f,GlobalId,'naming_authority','org.renci');
netcdf.putAtt(f,GlobalId,'id','ec2015');
netcdf.putAtt(f,GlobalId,'creator_name','Brian Blanton');
netcdf.putAtt(f,GlobalId,'creator_email','bblanton@renci.org');
netcdf.putAtt(f,GlobalId,'publisher_name','Brian Blanton');
netcdf.putAtt(f,GlobalId,'publisher_email','bblanton@renci.org');
netcdf.putAtt(f,GlobalId,'summary','Western North Atlantic, Caribbean and Gulf of Mexico Tidal Database, 2015.');

nn=fgs.nn;
ne=fgs.ne;
nf=size(D.UMAJOR,2);
maxc=7;
nface=size(fgs.e,2);

% def mode

DimEle=netcdf.defDim(f,'ele',ne);
DimNod=netcdf.defDim(f,'node',nn);
DimFqs=netcdf.defDim(f,'nfreqs',nf); 
DimChr=netcdf.defDim(f,'charlen',maxc);
DimFac=netcdf.defDim(f,'nface',nface); 

% variables
VarFqs = netcdf.defVar(f,'Frequencies','double',DimFqs);
    netcdf.putAtt(f,VarFqs,'long_name','Number of Frequencies/Constituents');
    netcdf.putAtt(f,VarFqs,'standard_name','Number of Frequencies/Constituents');
    netcdf.putAtt(f,VarFqs,'units','cycles per sec');

VarFqN = netcdf.defVar(f,'Names','char',[DimChr DimFqs]);
    netcdf.putAtt(f,VarFqN,'long_name','Constituent Names');
    
VarMaj = netcdf.defVar(f,'UMAJOR','double',[DimNod DimFqs]);
    netcdf.putAtt(f,VarMaj,'long_name','Tidal ellipse major axis length');
    netcdf.putAtt(f,VarMaj,'standard_name','Tidal ellipse major axis length');
    netcdf.putAtt(f,VarMaj,'units','m');
    netcdf.putAtt(f,VarMaj,'location','node');
    netcdf.putAtt(f,VarMaj,'coordinates','lon lat');
    netcdf.putAtt(f,VarMaj,'mesh','adcirc_mesh');

VarMin = netcdf.defVar(f,'UMINOR','double',[DimNod DimFqs]);
    netcdf.putAtt(f,VarMin,'long_name','Tidal ellipse minor axis length');
    netcdf.putAtt(f,VarMin,'standard_name','Tidal ellipse minor axis length');
    netcdf.putAtt(f,VarMin,'units','m');
    netcdf.putAtt(f,VarMin,'location','node');
    netcdf.putAtt(f,VarMin,'coordinates','lon lat');
    netcdf.putAtt(f,VarMin,'mesh','adcirc_mesh');
    
VarPha = netcdf.defVar(f,'PHASE','double',[DimNod DimFqs]);
    netcdf.putAtt(f,VarPha,'long_name','Tidal Ellipse phase [Time of Max Currents]');
    netcdf.putAtt(f,VarPha,'standard_name','Tidal Ellipse phase [Time of Max Currents]');
    netcdf.putAtt(f,VarPha,'units','radians');    
    netcdf.putAtt(f,VarPha,'location','node');
    netcdf.putAtt(f,VarPha,'coordinates','lon lat');
    netcdf.putAtt(f,VarPha,'mesh','adcirc_mesh');
    
VarOri = netcdf.defVar(f,'ORIEN','double',[DimNod DimFqs]);
    netcdf.putAtt(f,VarOri,'long_name','Tidal ellipse inclination');
    netcdf.putAtt(f,VarOri,'standard_name','Tidal ellipse inclination');
    netcdf.putAtt(f,VarOri,'units','radians CCW from east');
    netcdf.putAtt(f,VarOri,'location','node');
    netcdf.putAtt(f,VarOri,'coordinates','lon lat');
    netcdf.putAtt(f,VarOri,'mesh','adcirc_mesh');

VarLon = netcdf.defVar(f,'lon','double',DimNod);
    netcdf.putAtt(f,VarLon,'long_name','Longitude');
    netcdf.putAtt(f,VarLon,'standard_name','longitude');
    netcdf.putAtt(f,VarLon,'units','degrees_east');    
    netcdf.putAtt(f,VarLon,'location','node');
    netcdf.putAtt(f,VarLon,'axis','X');
    netcdf.putAtt(f,VarLon,'coordinates','lon lat');
    netcdf.putAtt(f,VarLon,'mesh','adcirc_mesh');
    
VarLat = netcdf.defVar(f,'lat','double',DimNod);
    netcdf.putAtt(f,VarLat,'long_name','Latitude');
    netcdf.putAtt(f,VarLat,'standard_name','latitude');
    netcdf.putAtt(f,VarLat,'units','degrees_north');    
    netcdf.putAtt(f,VarLat,'location','node');
    netcdf.putAtt(f,VarLat,'axis','Y');
    netcdf.putAtt(f,VarLat,'coordinates','lon lat');
    netcdf.putAtt(f,VarLat,'mesh','adcirc_mesh');
       
 VarDep = netcdf.defVar(f,'depth','double',DimNod);
    netcdf.putAtt(f,VarDep,'long_name','Bathymetry');
    netcdf.putAtt(f,VarDep,'standard_name','depth');
    netcdf.putAtt(f,VarDep,'units','m');    
    netcdf.putAtt(f,VarDep,'location','node');
    netcdf.putAtt(f,VarDep,'coordinates','lon lat');
    netcdf.putAtt(f,VarDep,'mesh','adcirc_mesh');
         
VarEle = netcdf.defVar(f,'ele','double',[DimEle DimFac]);
    netcdf.putAtt(f,VarEle,'long_name','element');
    netcdf.putAtt(f,VarEle,'start_index','1');
    netcdf.putAtt(f,VarEle,'cf_role','face_node_connectivity');
    
    
% VarRe = netcdf.defVar(f,'Re','double',[DimNod DimFqs]);
%     netcdf.putAtt(f,VarRe,'long_name','Real part of complex tide');
%     netcdf.putAtt(f,VarRe,'standard_name','Real part of complex tide');
%     netcdf.putAtt(f,VarRe,'units','m');
%  
%         
% VarIm = netcdf.defVar(f,'Im','double',[DimNod DimFqs]);
%     netcdf.putAtt(f,VarIm,'long_name','Imaginary part of complex tide');
%     netcdf.putAtt(f,VarIm,'standard_name','Imaginary part of complex tide');
%     netcdf.putAtt(f,VarIm,'units','m');
%     
netcdf.endDef(f);

names=D.PERNAMES;
for i=1:nf
    names{i}=strrep(names{i},'(','');
    names{i}=strrep(names{i},')','');
end

for i=1:nf
     l=length(names{i});
     netcdf.putVar(f,VarFqN,[0 i-1],[l 1],names{i});
end

% ReIm=F53.AMP.*exp(F53.PHA*pi/180);
% Re=real(ReIm);
% Im=imag(ReIm);
    
netcdf.putVar(f,VarFqs,0,nf,D.FREQ);
netcdf.putVar(f,VarMaj,[0 0],[nn nf],D.UMAJOR);
netcdf.putVar(f,VarMin,[0 0],[nn nf],D.UMINOR);
netcdf.putVar(f,VarPha,[0 0],[nn nf],D.PHASE);
netcdf.putVar(f,VarOri,[0 0],[nn nf],D.ORIEN);

netcdf.putVar(f,VarLon,0,nn,fgs.x);
netcdf.putVar(f,VarLat,0,nn,fgs.y);
netcdf.putVar(f,VarDep,0,nn,fgs.z);
netcdf.putVar(f,VarEle,[0 0],[ne nface],fgs.e);

% netcdf.putVar(f,VarRe,[0 0],[nn nf],Re);
% netcdf.putVar(f,VarIm,[0 0],[nn nf],Im);


netcdf.close(f);

    