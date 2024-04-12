function grs=LoadNWS13GroupNcFile(url)

if nargin==0
    error('Need url to group file')
end

% test file
% url='/Volumes/ees/MAPP/nws13test_Sept2018_test8.nc';
nc=ncgeodataset(url);

groupnames=split(nc.attribute{'group_order'});
ngroups=length(groupnames);
varnames={'lon' 'lat' 'PSFC' 'U10' 'V10' 'time'};
clear grs
grs=struct([]);

for i=1:ngroups
    grs(i).name=groupnames{i};
    for j=1:length(varnames)
        vstr=sprintf('%s/%s',groupnames{i},varnames{j});
        grs(i).(varnames{j})=nc{vstr};
        if strcmp(varnames{j},'time')
            ncc=nc{vstr};
            units=split(ncc.attribute('units'));
            stime=datetime(units{3});
            times=ncc.data(:);
            switch lower(units{1})
                case 'minutes'            
                    grs(i).time=stime+minutes(times);
                case 'seconds'            
                    grs(i).time=stime+seconds(times);
                otherwise
                    error(sprintf('uncoded time units: %s',units{1}))
            end
        end
    end
end

