function times=nctime(nc)
if any(ismember(nc.variables,'date'))
    times=datetime(datevec(nc.time('date')));
elseif any(ismember(nc.variables,'time'))
    times=datetime(datevec(nc.time('time')));
else
    error('Time variable not found in nc.')
end

