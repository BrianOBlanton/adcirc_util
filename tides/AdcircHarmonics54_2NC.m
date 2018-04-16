function AdcircHarmonics54_2NC(fgs,F54)

% append the grid like this: ncks -A ec2012_v3d.nc f54.nc

% output a netcdf file
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLOBBER'));
f = netcdf.create('f54.nc', mode);
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
nf=size(F54.FREQ,2);
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
    
VarUAmp = netcdf.defVar(f,'UAmp','double',[DimNod DimFqs]);
    netcdf.putAtt(f,VarUAmp,'long_name','EW Tide speed magnitude');
    netcdf.putAtt(f,VarUAmp,'standard_name','East/west tidal speed magnitude');
    netcdf.putAtt(f,VarUAmp,'units','m');
    netcdf.putAtt(f,VarUAmp,'location','node');
    netcdf.putAtt(f,VarUAmp,'coordinates','lon lat');
    netcdf.putAtt(f,VarUAmp,'mesh','adcirc_mesh');

VarUPha = netcdf.defVar(f,'UPha','double',[DimNod DimFqs]);
    netcdf.putAtt(f,VarUPha,'long_name','EW Tide speed phase');
    netcdf.putAtt(f,VarUPha,'standard_name','East/west tidal speed phase');
    netcdf.putAtt(f,VarUPha,'units','deg GMT');
    netcdf.putAtt(f,VarUPha,'location','node');
    netcdf.putAtt(f,VarUPha,'coordinates','lon lat');
    netcdf.putAtt(f,VarUPha,'mesh','adcirc_mesh');
    
VarVAmp = netcdf.defVar(f,'VAmp','double',[DimNod DimFqs]);
    netcdf.putAtt(f,VarVAmp,'long_name','NS Tide speed magnitude');
    netcdf.putAtt(f,VarVAmp,'standard_name','North/south tidal speed magnitude');
    netcdf.putAtt(f,VarVAmp,'units','m');
    netcdf.putAtt(f,VarVAmp,'location','node');
    netcdf.putAtt(f,VarVAmp,'coordinates','lon lat');
    netcdf.putAtt(f,VarVAmp,'mesh','adcirc_mesh');

VarVPha = netcdf.defVar(f,'VPha','double',[DimNod DimFqs]);
    netcdf.putAtt(f,VarVPha,'long_name','EW Tide speed phase');
    netcdf.putAtt(f,VarVPha,'standard_name','North/south tidal speed phase');
    netcdf.putAtt(f,VarVPha,'units','deg GMT');    
    netcdf.putAtt(f,VarVPha,'location','node');
    netcdf.putAtt(f,VarVPha,'coordinates','lon lat');
    netcdf.putAtt(f,VarVPha,'mesh','adcirc_mesh');
    
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

names=F54.PERNAMES;
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
    
netcdf.putVar(f,VarFqs,0,nf,F54.FREQ);
netcdf.putVar(f,VarUAmp,[0 0],[nn nf],F54.UA);
netcdf.putVar(f,VarUPha,[0 0],[nn nf],F54.UP);
netcdf.putVar(f,VarVAmp,[0 0],[nn nf],F54.VA);
netcdf.putVar(f,VarVPha,[0 0],[nn nf],F54.VP);

netcdf.putVar(f,VarLon,0,nn,fgs.x);
netcdf.putVar(f,VarLat,0,nn,fgs.y);
netcdf.putVar(f,VarDep,0,nn,fgs.z);
netcdf.putVar(f,VarEle,[0 0],[ne nface],fgs.e);

% netcdf.putVar(f,VarRe,[0 0],[nn nf],Re);
% netcdf.putVar(f,VarIm,[0 0],[nn nf],Im);


netcdf.close(f);

    