function TheGrid=ExtractGrid(NcTBHandle)
% Extract a grid from a CF-UGRID netCDF file
% TheGrid=ExtractGrid(NcTBHandle) 

verbose=false;

if ~isa(NcTBHandle,'ncgeodataset')
    error('Arg to ExtractGrid must be an ncgeodataset object from nctoolbox.')
end

TheGrid.name=NcTBHandle.attribute('agrid');
if isempty(TheGrid.name)
    TheGrid.name='unknown';
end

if any(strcmp(NcTBHandle.variables,'element'))
    TheGrid.e=double(NcTBHandle.data('element'));
elseif any(strcmp(NcTBHandle.variables,'ele'))
    TheGrid.e=double(NcTBHandle.data('ele'));
else
    error('Could not find an element list in variables.')
end

if size(TheGrid.e,1)<size(TheGrid.e,2)  % then element array is 3x, not x3
    TheGrid.e=TheGrid.e';
end

ConvertToCart=false;
temp=NcTBHandle.standard_name('longitude');
if ~isempty(temp)
    TheGrid.x=NcTBHandle.data(temp);
    ConvertToCart=true;
else
    if verbose, fprintf('*** No zonal variable with standard_name=longitude found.\n****** Looking for x_coordinate...\n'); end
    temp=NcTBHandle.standard_name('x_coordinate');
    if ~isempty(temp)
        if verbose, fprintf('****** Got it.\n');end
        TheGrid.x=NcTBHandle.data(temp);
    else
        error('\nNo zonal variable found with standard_name = {longitude,x_coordinate}')
    end
end

temp=NcTBHandle.standard_name('latitude');
if ~isempty(temp)
    TheGrid.y=NcTBHandle.data(temp);
else
    if verbose, fprintf('*** No meridional variable with standard_name=latitude found.\n****** Looking for y_coordinate...\n');end
    temp=NcTBHandle.standard_name('y_coordinate');
    if ~isempty(temp)
        if verbose, fprintf('****** Got it.\n');end
        TheGrid.y=NcTBHandle.data(temp);
    else
        error('\nNo meridional variable found with standard_name = {latitude,y_coordinate}')
    end
end

temp1=NcTBHandle.standard_name('depth_below_geoid');
temp2=NcTBHandle.standard_name('depth below geoid');
if ~(isempty(temp1) && isempty(temp2))
    temp=[temp1 temp2];
    temp=NcTBHandle.data(temp);
    TheGrid.z=cast(temp(:),'double');
else
    if verbose, fprintf('**** No depth variable with standard_name=depth_below_geoid found.  Setting depths to NaN...\n');end
    TheGrid.z=NaN(size(TheGrid.x));
end
       
if isa(TheGrid.x,'single')
    TheGrid.x=cast(TheGrid.x,'double');
    TheGrid.y=cast(TheGrid.y,'double');
    TheGrid.z=cast(TheGrid.z,'double');
end

TheGrid.bnd=detbndy(TheGrid.e);
TheGrid=el_areas(TheGrid);
TheGrid=attach_elem_centroids(TheGrid);
TheGrid=belint(TheGrid);
TheGrid.nn=size(TheGrid.x,1);
TheGrid.ne=size(TheGrid.e,1);

if ConvertToCart
    lon=TheGrid.x;
    lat=TheGrid.y;
    
    temp=TheGrid;
    
    [temp.x,temp.y]=AdcircCppForward(lon,lat,mean(lon),mean(lat));
    
    %fprintf('**** Lon/Lat grid converted to CPP. \n')
    temp=el_areas(temp);
    temp=belint(temp);
    temp=attach_elem_centroids(temp);
    
    TheGrid.ar_cart=temp.ar;
    TheGrid.A_cart=temp.A;
    TheGrid.A0_cart=temp.A0;
    TheGrid.B_cart=temp.B;
    TheGrid.T_cart=temp.T;
    TheGrid.dx_cart=temp.dx;
    TheGrid.dy_cart=temp.dy;
    TheGrid.dl_cart=sqrt(4*TheGrid.ar_cart/sqrt(3));
    TheGrid.x_cart=temp.x;
    TheGrid.y_cart=temp.y;
    
else
    
    TheGrid.ar_cart=TheGrid.ar;
    TheGrid.A_cart=TheGrid.A;
    TheGrid.A0_cart=TheGrid.A0;
    TheGrid.B_cart=TheGrid.B;
    TheGrid.T_cart=TheGrid.T;
    TheGrid.dx_cart=TheGrid.dx;
    TheGrid.dy_cart=TheGrid.dy;
    TheGrid.dl_cart=sqrt(4*TheGrid.ar_cart/sqrt(3));
    TheGrid.x_cart=TheGrid.x;
    TheGrid.y_cart=TheGrid.y;
    
end
    
%TheGrid=el_areas(TheGrid);
if all(TheGrid.ar<0)  % assume elements are ordered CW and switch to CCW
    fprintf('**** Permuting element list to CCW.\n ')
    TheGrid.e=TheGrid.e(:,[1 3 2]);
    TheGrid=el_areas(TheGrid);
    TheGrid.bnd=detbndy(TheGrid.e);
end


%
%LabSig  Brian O. Blanton
%        Renaissance Computing Institute
%        University of North Carolina
%        Chapel Hill, NC
%
%        brian_blanton@renci.org
%

