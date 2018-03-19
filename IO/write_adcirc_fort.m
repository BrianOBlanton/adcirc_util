function write_adcirc_fort(fortdatastruct,filename)
% call as:  write_adcirc_fort(fortdatastruct,filename);

% if exist(filename,'file')
%     error('%s already exists. Delete and retry.',filename)
% end


if isfield(fortdatastruct,'zeta')
    [m,n]=size(fortdatastruct.zeta);
    data=fortdatastruct.zeta;
%    data(isnan(data))=-99999;
    nf=1;
elseif isfield(fortdatastruct,'ubar')
    [m,n]=size(fortdatastruct.ubar);
    datau=fortdatastruct.ubar;
    datav=fortdatastruct.vbar;
    datau(isnan(datau))=-99999;
    datav(isnan(datav))=-99999;
    nf=2;
else
    error('fortdatastruct does not contain a needed field');
end


fid=fopen(filename,'w');

% niter=(fortdatastruct.t(end)-fortdatastruct.t(1))/(n-1);
% if n==1, niter=1;,end

fprintf(fid,'%s\n',deblank(fortdatastruct.header));
fprintf(fid,'%10d %10d %f %10d %2d FileFmtVersion:    1050624\n',n,m,fortdatastruct.dt,fortdatastruct.NSTEP,fortdatastruct.IFLAG);

for i=1:n
    fprintf(fid,'       %16.10e %12d\n',fortdatastruct.t(i),fortdatastruct.iter(i));
    if nf==1
        out=[(1:m)' data(:,i)];
        fprintf(fid,'%12d    %16.10e\n',out');
    else
        out=[(1:m)' datau(:,i)  datav(:,i)];
        fprintf(fid,'%12d    %16.10e   %16.10e\n',out');
    end
    
end


fclose(fid);
