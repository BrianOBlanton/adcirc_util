function TheGrid=ExtractGrid(NcTBHandle)
% TheGrid=ExtractGrid(NcTBHandle)

if ~isa(NcTBHandle,'ncgeodataset')
    error('Arg to ExtractGrid must be an ncgeodataset object from nctoolbox.')
end

TheGrid.name=NcTBHandle.attribute('agrid');
if isempty(TheGrid.name)
    TheGrid.name='unknown';
end

TheGrid.e=double(NcTBHandle.data('element'));

temp=NcTBHandle.standard_name('longitude');
if ~isempty(temp)
    TheGrid.x=NcTBHandle.data(temp);
else
    error('**** No x-coord variable with standard_name=longitude found.')
end

temp=NcTBHandle.standard_name('latitude');
if ~isempty(temp)
    TheGrid.y=NcTBHandle.data(temp);     
else
    error('**** No y-coord variable with standard_name=latitude found.')
end

temp1=NcTBHandle.standard_name('depth_below_geoid');
temp2=NcTBHandle.standard_name('depth below geoid');
if ~isempty(temp1 || temp2)
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
