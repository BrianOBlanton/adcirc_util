function err=Owi2NetCdf(FileName)
% Owi2NetCdf writes OWI win or pre file into netCDF for nws14.  
% Can't (yet) write a single nc file with both wind and pressure.
% If merged is needed, use nco's ncks
% err=Owi2NetCdf(<FileName>.221|222|223|224|win|wnd|pre)

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
            ncwriteatt(f_out,'time','calendar','gregorian');
            %ncwriteatt(f_out,'time','start_index',1);

        nccreate(f_out,'longitude','Dimensions',{'longitude', iLong},'Format','netcdf4');
            ncwriteatt(f_out,'longitude','long_name','longitude');
            ncwriteatt(f_out,'longitude','standard_name','longitude');
            ncwriteatt(f_out,'longitude','units','degrees_east');
            ncwriteatt(f_out,'longitude','positive','east');
            ncwriteatt(f_out,'longitude','axis','X');

        nccreate(f_out,'latitude','Dimensions',{'latitude', iLat},'Format','netcdf4');
            ncwriteatt(f_out,'latitude','long_name','latitude');
            ncwriteatt(f_out,'latitude','standard_name','latitude');
            ncwriteatt(f_out,'latitude','units','degrees_north');
            ncwriteatt(f_out,'latitude','positive','north');
            ncwriteatt(f_out,'latitude','axis','Y');

        if any(strcmpi(Type,{'222','224','wnd','win'})) 
        
            nccreate(f_out,'u10','Dimensions',{'longitude' 'latitude' 'time'},'Format','netcdf4','ChunkSize',[iLong iLat 1]);
                ncwriteatt(f_out,'u10','long_name','10 meter U wind component');
                ncwriteatt(f_out,'u10','standard_name','eastward_wind');
                ncwriteatt(f_out,'u10','units','m/s');
                ncwriteatt(f_out,'u10','positive','east');
                ncwriteatt(f_out,'u10','coordinates','longitude latitude time');

            nccreate(f_out,'v10','Dimensions',{'longitude' 'latitude' 'time'},'Format','netcdf4','ChunkSize',[iLong iLat 1]);
                ncwriteatt(f_out,'v10','long_name','10 meter V wind component');
                ncwriteatt(f_out,'v10','standard_name','northward_wind');
                ncwriteatt(f_out,'v10','units','m/s');
                ncwriteatt(f_out,'v10','positive','north');
                ncwriteatt(f_out,'v10','coordinates','longitude latitude time');
            
        else
           nccreate(f_out,'msl','Dimensions',{'longitude' 'latitude' 'time'},'Format','netcdf4','ChunkSize',[iLong iLat 1]);
                ncwriteatt(f_out,'msl','long_name','Mean sea level pressure');
                ncwriteatt(f_out,'msl','standard_name','air_pressure_at_mean_sea_level');
                ncwriteatt(f_out,'msl','units','Pa');
                ncwriteatt(f_out,'msl','positive','east');
                ncwriteatt(f_out,'msl','coordinates','longitude latitude time');
        end
           
        x=A(6)+(0:iLong-1)*A(3);
        y=A(5)+(0:iLat-1)*A(4);
        ncwrite(f_out,'longitude', x); 
        ncwrite(f_out,'latitude', y);
        
        starttime=ctime;
        
    end
   
    if Debug
        if any(strcmpi(Type,{'222','224','wnd','win'})) 
            t='u10';
            fprintf('Scanning [%d x %d] %s values at time=%s\n',iLong,iLat,t,datestr(ctime,0))
        end
    end
    u=fscanf(fid,'%f',[iLong,iLat]);  

    if any(strcmpi(Type,{'222','224','wnd','win'}))
        if Debug
            fprintf('Scanning [%d x %d] v10 values at time=%s\n',iLong,iLat,datestr(ctime,0))
        end
        v=fscanf(fid,'%f',[iLong,iLat]);
    end

    t= round((ctime-starttime)*86400);    
    if Debug, fprintf('Flushing data to nc file at t=%d ... ',t);end

    ncwrite(f_out,'time',  t, k);
   
    if any(strcmpi(Type,{'222','224','wnd','win'})) 
        ncwrite(f_out,'u10', u, [1 1 k]); 
        ncwrite(f_out,'v10', v, [1 1 k]); 
    else
        % convert mb to Pa
        ncwrite(f_out,'msl', u*100, [1 1 k]); 
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