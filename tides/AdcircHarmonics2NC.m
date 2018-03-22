function AdcircHarmonics2NC(fgs,F53)



% output a netcdf file
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLOBBER'));
f = netcdf.create('f53.nc', mode);
GlobalId=netcdf.getConstant('GLOBAL');
netcdf.putAtt(f,GlobalId,'creation_date',datestr(now));
netcdf.putAtt(f,GlobalId,'title','Tide Harmonics');
netcdf.putAtt(f,GlobalId,'model_domain',strrep(fgs.name,'.grd',''));



nn=fgs.nn;
nf=size(F53.FREQ,2);
maxc=7;

% def mode

DimNod=netcdf.defDim(f,'node',nn);
DimFqs=netcdf.defDim(f,'nfreqs',nf); 
DimChr=netcdf.defDim(f,'charlen',maxc);

% variables
VarFqs = netcdf.defVar(f,'Frequencies','double',DimFqs);
    netcdf.putAtt(f,VarFqs,'long_name','Number of Frequencies/Constituents');
    netcdf.putAtt(f,VarFqs,'standard_name','Number of Frequencies/Constituents');
    netcdf.putAtt(f,VarFqs,'units','cycles per sec');

VarFqN = netcdf.defVar(f,'Names','char',[DimChr DimFqs]);
    netcdf.putAtt(f,VarFqN,'long_name','Constituent Names');
    
VarAmp = netcdf.defVar(f,'Amp','double',[DimNod DimFqs]);
    netcdf.putAtt(f,VarAmp,'long_name','Tide Amplitudes');
    netcdf.putAtt(f,VarAmp,'standard_name','Tide Amplitudes');
    netcdf.putAtt(f,VarAmp,'units','m');

VarPha = netcdf.defVar(f,'Pha','double',[DimNod DimFqs]);
    netcdf.putAtt(f,VarPha,'long_name','Tide Phases');
    netcdf.putAtt(f,VarPha,'standard_name','Tide Phases');
    netcdf.putAtt(f,VarPha,'units','deg GMT');
    
    
VarRe = netcdf.defVar(f,'Re','double',[DimNod DimFqs]);
    netcdf.putAtt(f,VarRe,'long_name','Real part of complex tide');
    netcdf.putAtt(f,VarRe,'standard_name','Real part of complex tide');
    netcdf.putAtt(f,VarRe,'units','m');
 
        
VarIm = netcdf.defVar(f,'Im','double',[DimNod DimFqs]);
    netcdf.putAtt(f,VarIm,'long_name','Imaginary part of complex tide');
    netcdf.putAtt(f,VarIm,'standard_name','Imaginary part of complex tide');
    netcdf.putAtt(f,VarIm,'units','m');
    
netcdf.endDef(f);

%netcdf.putVar(f,VarFqN,[0 0],[maxc nf],char(F53.PERNAMES(:))');

for i=1:nf
     l=length(F53.PERNAMES{i});
     netcdf.putVar(f,VarFqN,[0 i-1],[l 1],F53.PERNAMES{i});
end


ReIm=F53.AMP.*exp(F53.PHA*pi/180);
Re=real(ReIm);
Im=imag(ReIm);
    
netcdf.putVar(f,VarFqs,0,nf,F53.FREQ);
netcdf.putVar(f,VarAmp,[0 0],[nn nf],F53.AMP);
netcdf.putVar(f,VarPha,[0 0],[nn nf],F53.PHA);
netcdf.putVar(f,VarRe,[0 0],[nn nf],Re);
netcdf.putVar(f,VarIm,[0 0],[nn nf],Im);


netcdf.close(f);

    