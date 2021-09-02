function D=read_adcirc_fort62(fname,asd)
%READ_ADCIRC_FORT62 read ADCIRC station velocity output file
% READ_ADCIRC_FORT62 reads in the contents of an ADCIRC station 
% velocity file, typically output to IO unit 62 (fort.62).  
% This output file is described in the ADCIRC Users Manual as 
% "Depth-averaged Velocity Time Series at Specified Velocity 
% Recording Stations (fort.62)". Time is converted from ADCIRC time 
% to Gregorian time.
%
%  Input : fname - filename from which to read elevation.  If empty,
%                  routine reads from fort.62. 
%          asd   - ADCIRC start date, needed for conversion of 
%                  ADCIRC time to Gregorian (quoddy) time.  If empty,
%                  no conversion is made.  ASD must be passed in as a
%                  1x6-element vector, [yyyy,mm,dd,hr,mn,sec].
%
% Output : a struct array is returned with the following fields:
%          D.u     - matrix of E/W velocity, one station per column
%          D.v     - matrix of N/S velocity, one station per column
%          D.t     - vector of times
%          D.nstae - number of stations
%          D.dt    - output interval in fort.62 file
%
% Call as: D=read_adcirc_fort62(fname,asd);
% 
% Written by: Brian Blanton, Spring '03


if nargin==0 & nargout==0
  disp('D=read_adcirc_fort62(fname,asd);')
  return
else
   if nargin==1
   % See if fname is string
      if isempty(fname)
         fname='fort.62';
      elseif ~isstr(fname)
         error('Filename to READ_ADCIRC_FORT62 must be a string.')
      end
      time_offset=0.;
   elseif nargin==2
   % See if fname is string
      if isempty(fname)
         fname='fort.62';
      elseif ~isstr(fname)
         error('Filename to READ_ADCIRC_FORT62 must be a string.')
      end
      if isempty(asd)
         time_offset=0.;      
      else
         [m,n]=size(asd);
         if n~=6
            error('ASD is not correctly sized.')
         else
            time_offset=datenum(asd);
         end
      end
    else
       time_offset=0.;      
       fname='fort.62';
   end
end

disp(['Gregorian Time Offset = ' num2str(time_offset)])

[fid,message]=fopen(fname,'r');
if fid==-1
   error(['Could not open ' fname ' because ' message])
end

% The header line: RUNDES, RUNID, AGRID 
header=fgets(fid);

% NTRSPE, NSTAE, DT*NSPOOLE, NSPOOLE, IRTYPE 
data=fscanf(fid,'%d %d %f %d %d',[5])';
fgets(fid);

NTRSPE=data(1);
NSTAE=data(2);
outdt=data(3);
NSPOOLE=data(4);

u=NaN*ones(NTRSPE,NSTAE);
v=NaN*ones(NTRSPE,NSTAE);
t=NaN*ones(NTRSPE,1);

for ii=1:NTRSPE
   if feof(fid),break,end
   temp=fscanf(fid,'%f %d',[1 2]);
   t(ii)=temp(1);
   data=fscanf(fid,'%d %f %f',[3 NSTAE]);
   u(ii,:)=data(2,:);
   v(ii,:)=data(3,:);
   fgets(fid);
end

D.u=u;
D.v=v;
D.t=t;
if time_offset>0
   D.t=D.t/86400+time_offset;
end
D.nstae=NSTAE;
D.dt=outdt;

%
%LabSig  Brian O. Blanton
%        Department of Marine Sciences
%        12-7 Venable Hall
%        CB# 3300
%        University of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        brian_blanton@unc.edu
%
