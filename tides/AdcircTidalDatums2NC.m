function AdcircTidalDatums2NC(fgs,D)


% output a netcdf file
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLOBBER'));
f = netcdf.create('tideDatums.nc', mode);
GlobalId=netcdf.getConstant('GLOBAL');
netcdf.putAtt(f,GlobalId,'creation_date',datestr(now));
netcdf.putAtt(f,GlobalId,'title','Tide Datums');
netcdf.putAtt(f,GlobalId,'model_domain',strrep(fgs.name,'.grd',''));

nn=fgs.nn;

% def mode

DimNod=netcdf.defDim(f,'node',nn);
DimDTS=netcdf.defDim(f,'ndats',4); 


% variables
VarMhhw = netcdf.defVar(f,'MHHW','double',DimNod);
    netcdf.putAtt(f,VarMhhw,'long_name','Mean Higher High Water');
    netcdf.putAtt(f,VarMhhw,'units','m');
    netcdf.putAtt(f,VarMhhw,'standard_name','MHHW');

VarMllw = netcdf.defVar(f,'MLLW','double',DimNod);
    netcdf.putAtt(f,VarMllw,'long_name','Mean Lower Low Water');
    netcdf.putAtt(f,VarMllw,'units','m');
    netcdf.putAtt(f,VarMllw,'standard_name','MLLW');

VarHat = netcdf.defVar(f,'HAT','double',DimNod);
    netcdf.putAtt(f,VarHat,'long_name','Highest Astronomical Tide');
    netcdf.putAtt(f,VarHat,'units','m');
    netcdf.putAtt(f,VarHat,'standard_name','HAT');
   
VarSpT = netcdf.defVar(f,'SpringTide','double',DimNod);
    netcdf.putAtt(f,VarSpT,'long_name','Spring Tide');
    netcdf.putAtt(f,VarSpT,'units','m');
    netcdf.putAtt(f,VarSpT,'standard_name','Spring Tide');
    
VarNpT = netcdf.defVar(f,'NeapTide','double',DimNod);
    netcdf.putAtt(f,VarNpT,'long_name','Neap Tide');
    netcdf.putAtt(f,VarNpT,'units','m');
    netcdf.putAtt(f,VarNpT,'standard_name','Neap Tide');
    
netcdf.endDef(f);

netcdf.putVar(f,VarMhhw,0,nn,D.MHHW);
netcdf.putVar(f,VarMllw,0,nn,D.MLLW);
netcdf.putVar(f,VarMllw,0,nn,D.HAT);
netcdf.putVar(f,VarSpT,0,nn,D.SpringTide);
netcdf.putVar(f,VarNpT,0,nn,D.NeapTide);


netcdf.close(f);

