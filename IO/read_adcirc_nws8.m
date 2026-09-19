function D=read_adcirc_nws8(filename)
% read_adcirc_nws8 - Read an ADCIRC NWS8 fort.22 file
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

%BN, CY, YYYYMMDDHH, MIN,TECH, TAU,latNS, lonEW,VMAX, MSLP, TY, RAD,  WC, RAD1, RAD2, RAD3, RAD4, POUT, ROUT, RMW,GUST, EYE,  SR,MAXS,INIT, DIR, SPD,  STORMNAME

[D.basin,D.stormnumber,yyyymmdd,blnk,blnk,hr,lat,lon,rmw,cp]=textread(filename,'%s%d%s%s%s%d%s%s%s%s%*[^\n]','delimiter',',');


for i=1:length(lat)
   D.lat(i)=str2double(replace(lat{i},{'N','S'},''))/10;
   D.lon(i)=-str2double(replace(lon{i},{'E','W'},''))/10;
   D.rmw(i)=str2double(rmw{i});
   D.cp(i)=str2double(cp{i});
   yyyy=str2num(yyyymmdd{i}(1:4));
   mm=str2num(yyyymmdd{i}(5:6));
   dd=str2num(yyyymmdd{i}(7:8));
   D.time(i)=datenum(yyyy,mm,dd,hr(i),0,0);
end

