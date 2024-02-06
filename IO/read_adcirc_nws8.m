function D=read_adcirc_nws9(filename)
% read_adcirc_nws9 - Read an ADCIRC NWS8 fort.22 file
% returns a struct with: 
% 
% f22 = 
%          basin: 
%    stormnumber:
%            lat:
%            lon:
%            rmw:
%             cp:
%           time:

[D.basin,D.stormnumber,yyyymmdd,blnk,blnk,hr,lat,lon,rmw,cp]=textread(filename,'%s%d%s%s%s%d%s%s%s%s%*[^\n]','delimiter',',');


for i=1:length(lat)
   D.lat(i)=str2num(lat{i}(1:3))/10;
   D.lon(i)=-str2num(lon{i}(1:3))/10;
   D.rmw(i)=str2num(rmw{i});
   D.cp(i)=str2num(cp{i});
   yyyy=str2num(yyyymmdd{i}(1:4));
   mm=str2num(yyyymmdd{i}(5:6));
   dd=str2num(yyyymmdd{i}(7:8));
   D.time(i)=datenum(yyyy,mm,dd,hr(i),0,0);
end

