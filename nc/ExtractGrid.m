function TheGrid=ExtractGrid(NcTBHandle)
% Extract a grid from a CF-UGRID netCDF file
% TheGrid=ExtractGrid(NcTBHandle) 

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

NeedsConvertToCart=true;
temp=NcTBHandle.standard_name('longitude');
if ~isempty(temp)
    TheGrid.x=NcTBHandle.data(temp);
else
    fprintf('**** No zonal variable with standard_name=longitude found. Looking for x_coordinate...')
    temp=NcTBHandle.standard_name('x_coordinate');
    if ~isempty(temp)
        fprintf(' Got it.\n')
        TheGrid.x=NcTBHandle.data(temp);
    else
        error('\nNo zonal variable found with standard_name = {longitude,x_coordinate}')
    end
    NeedsConvertToCart=false;

end

temp=NcTBHandle.standard_name('latitude');
if ~isempty(temp)
    TheGrid.y=NcTBHandle.data(temp);
else
    fprintf('**** No meridional variable with standard_name=latitude found. Looking for y_coordinate...')
    temp=NcTBHandle.standard_name('y_coordinate');
    if ~isempty(temp)
        fprintf(' Got it.\n')
        TheGrid.y=NcTBHandle.data(temp);
    else
        error('\nNo meridional variable found with standard_name = {latitude,y_coordinate}')
    end
    NeedsConvertToCart=false;

end

if NeedsConvertToCart
    [TheGrid.x,TheGrid.y]=AdcircCppForward(TheGrid.x,TheGrid.y,-80,36);
    fprintf('**** Lon/Lat grid converted to CPP. ...')

end
    
temp1=NcTBHandle.standard_name('depth_below_geoid');
temp2=NcTBHandle.standard_name('depth below geoid');
if ~(isempty(temp1) && isempty(temp2))
    temp=[temp1 temp2];
    temp=NcTBHandle.data(temp);
    TheGrid.z=cast(temp(:),'double');
else
    sprintf('**** No depth variable with standard_name=depth_below_geoid found.  Setting depths to NaN...')
    TheGrid.z=NaN(size(TheGrid.x));
end
       
if isa(TheGrid.x,'single')
    TheGrid.x=cast(TheGrid.x,'double');
    TheGrid.y=cast(TheGrid.y,'double');
    TheGrid.z=cast(TheGrid.z,'double');
end

TheGrid.bnd=detbndy(TheGrid.e);

TheGrid=el_areas(TheGrid);
if all(TheGrid.ar<0)  % assume elements are ordered CW and switch to CCW
    SetUIStatusMessage('**** Permuting element list to CCW. ')
    TheGrid.e=TheGrid.e(:,[1 3 2]);
    TheGrid=el_areas(TheGrid);
    TheGrid.bnd=detbndy(TheGrid.e);
end
TheGrid=belint(TheGrid);

TheGrid=attach_elem_centroids(TheGrid);