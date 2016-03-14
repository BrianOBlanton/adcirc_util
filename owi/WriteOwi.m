function WriteOwi(D,filename)
%WriteOwi
% WriteOwi(D,filename)


vectorfile=0;
if isfield(D,'v')
   vectorfile=1;
end

if ~exist('filename')
   filename='test.owi.221';
   if vectorfile
      filename='test.owi.222';
   end
end

time_string='iLat=%4diLong=%4dDX=%6.4fDY=%6.4fSWLat=%8.5fSWLon=%8.4fDT=%12s';
%value_string=' %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n';
value_string=' %9d %9d %9d %9d %9d %9d %9d %9d\n';

fid=fopen(filename,'w');

t1=datestr(D(1).time,30);
t1([9 14 15])=[];
t2=datestr(D(end).time,30);
t2([9 14 15])=[];
header=sprintf('Oceanweather WIN/PRE Format                        %12s     %12s',t1,t2);
fprintf(fid,'%s\n',header);

for i=1:length(D)
   t=datestr(D(i).time,30);
   t([9 14 15])=[];
   header=sprintf(time_string,D(i).iLat,D(i).iLong,D(i).DX,D(i).DY,D(i).SWLat,D(i).SWLon,t);
   disp(sprintf('%s',header))
   fprintf(fid,'%s\n',header);
   out=D(i).u;
   disp(sprintf('   Min,Max u = %f %f',min(out(:)),max(out(:))))
   fprintf(fid,value_string,out);
   fprintf(fid,'\n');
   if vectorfile
      out=D(i).v;
      fprintf(fid,value_string,out);
      fprintf(fid,'\n');
      disp(sprintf('   Min,Max v = %f %f',min(out(:)),max(out(:))))
   end
   
end

fclose(fid);
