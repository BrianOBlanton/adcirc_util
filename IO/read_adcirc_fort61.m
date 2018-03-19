function D=read_adcirc_fort61(fname,asd)
%READ_ADCIRC_FORT61 read ADCIRC station elevation output file
% READ_ADCIRC_FORT61 reads in the contents of an ADCIRC station 
% elevation file, typically output to IO unit 61 (fort.61).  
% This output file is described in the ADCIRC Users Manual as 
% "Elevation Time Series at Specified Elevation Recording Stations (fort.61)". 
% Time is converted from ADCIRC time to Gregorian time.
%
%  Input : fname - filename from which to read elevation.  If empty,
%                  routine reads from fort.61. 
%          asd   - ADCIRC start date, needed for conversion of 
%                  ADCIRC time to Gregorian (quoddy) time.  If empty,
%                  no conversion is made.  ASD must be passed in as a
%                  1x6-element vector, [yyyy,mm,dd,hr,mn,sec].
%
% Output : a struct array is returned with the following fields:
%          D.zeta  - matrix of elevation, one station per column
%          D.t     - vector of times
%          D.nstae - number of stations
%          D.dt    - output interval in fort.61 file
%
% Call as: D=read_adcirc_fort61(fname,asd);
% 
% Written by: Brian Blanton, Spring '02


Clip=1;

if nargin==0 & nargout==0
  disp('D=read_adcirc_fort61(fname,asd);')
  return
else
   if nargin==1
   % See if fname is string
      if isempty(fname)
         fname='fort.61';
      elseif ~isstr(fname)
         error('Filename to READ_ADCIRC_FORT61 must be a string.')
      end
      time_offset=0.;
   elseif nargin==2
   % See if fname is string
      if isempty(fname)
         fname='fort.61';
      elseif ~isstr(fname)
         error('Filename to READ_ADCIRC_FORT61 must be a string.')
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
       fname='fort.61';
   end
end

%disp(['Gregorian Time Offset = ' num2str(time_offset)])

[fid,message]=fopen(fname,'r');
if fid==-1
   error(['Could not open ' fname ' because ' message])
end

% The header line: RUNDES, RUNID, AGRID 
header=fgets(fid);

% NTRSPE, NSTAE, DT*NSPOOLE, NSPOOLE, IRTYPE 
data=fscanf(fid,'%d %d %f %d %d',[5])';

NTRSPE=data(1);
NSTAE=data(2);
outdt=data(3);
NSPOOLE=data(4);

fgetl(fid);

zeta=NaN*ones(NTRSPE,NSTAE);
t=NaN*ones(NTRSPE,1);

for ii=1:NTRSPE
   if feof(fid),break,end
   temp=fscanf(fid,'%f %d',[1 2]);
   t(ii)=temp(1);
   data=fscanf(fid,'%d %f',[2 NSTAE]);
   zeta(ii,:)=data(2,:);
   fgets(fid);
end
fclose(fid);

D.zeta=zeta;
D.t=t;
if time_offset>0
   D.t=D.t/86400+time_offset;
end
D.nstae=NSTAE;
D.dt=outdt;

if Clip
   inan=find(D.zeta<-9999);
   D.zeta(inan)=NaN;
end


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
