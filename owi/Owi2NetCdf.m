function err=Owi2NetCdf(FileName)
% err=Owi2NetCdf(FileName)

err=0;

Debug=true;

if ~exist(FileName,'file')
   error('%s does not exist.',FileName)
end

f_out=[FileName '.nc'];
if exist(f_out,'file')
    error([f_out ' already exists. Terminal.'])
end

[~,~,Type]=fileparts(FileName);
Type=Type(2:end);
if ~any(strcmpi(Type,{'221','222','223','224','pre','wnd','win'}))
   error('WND/PRE file extension (%s) must be 221|222|223|224|win|wnd|pre.',Type)
end
fid=fopen(FileName,'r');

headerline=fgets(fid);

%             A(1)    A(2)  A(3) A(4)  A(5)    A(6)   A(7)
time_string='iLat=%diLong=%dDX=%fDY=%fSWLat=%fSWLon=%fDT=%4d%2d%2d%2d%2d';

k=0;

if Debug>-1
   fprintf('Reading from %s \n',FileName)
end

while ~feof(fid)
   k=k+1;
   l=fgets(fid);
   if Debug, fprintf('%d : %s',k,l);end
   %[iLat,iLong,DX,DY,SWLat,SWLon,y,m,d,h,mn]
   A=sscanf(l,time_string);
   ctime=datenum(A(7),A(8),A(9),A(10),A(11),0);
   
   iLat=A(1);
   iLong=A(2);

   if k==1
        fprintf('Creating %s\n',f_out);
        nccreate(f_out,'time','Dimensions',{'time', Inf},'Format','netcdf4');
            ncwriteatt(f_out,'time','long_name','model time');
            ncwriteatt(f_out,'time','standard_name','time');
            ncwriteatt(f_out,'time','units',sprintf('seconds since %s',datestr(ctime,'yyyy-mm-dd HH:MM:SS')));
            ncwriteatt(f_out,'time','start_index',1);

        nccreate(f_out,'x','Dimensions',{'x', iLong},'Format','netcdf4');
            ncwriteatt(f_out,'x','long_name','longitude');
            ncwriteatt(f_out,'x','standard_name','longitude');
            ncwriteatt(f_out,'x','units','degrees_east');
            ncwriteatt(f_out,'x','positive','east');

        nccreate(f_out,'y','Dimensions',{'y', iLat},'Format','netcdf4');
            ncwriteatt(f_out,'y','long_name','latitude');
            ncwriteatt(f_out,'y','standard_name','latitude');
            ncwriteatt(f_out,'y','units','degrees_north');
            ncwriteatt(f_out,'y','positive','north');

        if any(strcmpi(Type,{'222','224','wnd','win'})) 
        
            nccreate(f_out,'windx','Dimensions',{'y' 'x' 'time'},'Format','netcdf4'); % , 'ChunkSize',  [iLong iLat 1]);
                ncwriteatt(f_out,'windx','long_name','e/w wind velocity');
                ncwriteatt(f_out,'windx','standard_name','eastward_wind');
                ncwriteatt(f_out,'windx','units','m/s');
                ncwriteatt(f_out,'windx','positive','east');

            nccreate(f_out,'windy','Dimensions',{'y' 'x' 'time'},'Format','netcdf4'); % , 'ChunkSize',  [iLong iLat 1]);
                ncwriteatt(f_out,'windy','long_name','n/s wind velocity');
                ncwriteatt(f_out,'windy','standard_name','northward_wind');
                ncwriteatt(f_out,'windy','units','m/s');
                ncwriteatt(f_out,'windy','positive','north');
            
        else
           nccreate(f_out,'pressure','Dimensions',{'y' 'x' 'time'},'Format','netcdf4'); % , 'ChunkSize',  [iLong iLat 1]);
                ncwriteatt(f_out,'pressure','long_name','air pressure at sea level');
                ncwriteatt(f_out,'pressure','standard_name','air_pressure_at_sea_level');
                ncwriteatt(f_out,'pressure','units','meters of water');
                ncwriteatt(f_out,'pressure','positive','east');
        end
           
        x=A(6)+(0:iLong-1)*A(3);
        y=A(5)+(0:iLat-1)*A(4);
        ncwrite(f_out,'x', x); 
        ncwrite(f_out,'y', y);
        
        starttime=ctime;
        
    end
   
    if Debug
        if any(strcmpi(Type,{'222','224','wnd','win'})) 
            t='windx';
            fprintf('Scanning [%d x %d] %s values at time=%s\n',iLong,iLat,t,datestr(ctime,0))
        end
    end
    u=fscanf(fid,'%f',[iLong,iLat]);  

    if any(strcmpi(Type,{'222','224','wnd','win'}))
        if Debug
            fprintf('Scanning [%d x %d] windy values at time=%s\n',iLong,iLat,datestr(ctime,0))
        end
        v=fscanf(fid,'%f',[iLong,iLat]);
    end

    t= round((ctime-starttime)*86400);    
    if Debug, fprintf('Flushing data to nc file at t=%d ... ',t);end

    ncwrite(f_out,'time',  t, k);
    
    if any(strcmpi(Type,{'222','224','wnd','win'})) 
        ncwrite(f_out,'windx', u, [1 1 k]); 
        ncwrite(f_out,'windy', v, [1 1 k]); 
    else
        ncwrite(f_out,'pressure', u, [1 1 k]); 
    end
    if Debug, fprintf('done.\n');end

    if Debug
        fprintf('   Min,Max u = %f %f\n',min(u(:)),max(u(:)))
    end
    if any(strcmpi(Type,{'222','224','wnd','win'}))
        if Debug
            fprintf('   Min,Max v = %f %f\n',min(v(:)),max(v(:)))
        end
    end
    fgetl(fid);
    if Debug
        fprintf('\n'); 
    end
end

fclose(fid);


   
   

%% nc test
% delete myncclassic.nc
% nccreate('myncclassic.nc','peaks','Dimensions',{'r' 200 'c' 200 't' Inf},'Format','netcdf4');
% ncwrite('myncclassic.nc','peaks', peaks(200),[1 1 1]);
% ncdisp('myncclassic.nc');
         


