function D=read_adcirc_fort63(varargin)
%READ_ADCIRC_FORT63 read ADCIRC global elevation output file
% READ_ADCIRC_FORT63 reads in the contents of an ADCIRC global 
% elevation file, typically output to IO unit 63 (fort.63).  
% This output file is described in the ADCIRC Users Manual as 
% "Elevation Time Series at All Nodes in the Model Grid (fort.63)". 
% Time is converted from ADCIRC time to Gregorian time.
%
%  Input : fname - filename from which to read elevation.  If empty,
%                  routine reads from fort.63. 
%          asd   - ADCIRC start date, needed for conversion of 
%                  ADCIRC time to Gregorian (quoddy) time.  If empty,
%                  no conversion is made.  ASD must be passed in as a
%                  1x6-element vector, [yyyy,mm,dd,hr,mn,sec].
%
% Output : a struct array is returned with the following fields:
%          D.zeta  - matrix of elevation, one node per row
%          D.t     - vector of times
%          D.nodes - number of nodes
%          D.dt    - output interval in fort.63 file
%
% Call as: D=read_adcirc_fort63(FortName,<FortName>,...);
% 
% Written by: Brian Blanton, Spring '03
%             changed to varargins, Fall '07


if nargin==0 & nargout==0
   disp('D=read_adcirc_fort63(''FortName'',''<FortName>'',...);')
   return
end


% defaults
Compact=1;
ClipVal=NaN;
FortName='fort.63';
ASD=0;   % Adcirc Start Date
TimeClipMin=NaN;
TimeClipMax=NaN;


% Strip off propertyname/value pairs in varargin not related to
% "line" object properties.
k=1;
while k<length(varargin),
  switch lower(varargin{k}),
    case 'compact',
      Compact=varargin{k+1};
      varargin([k k+1])=[];
    case 'clipval',
      ClipVal=varargin{k+1};
      varargin([k k+1])=[];
    case 'fortname',
      FortName=varargin{k+1};
      varargin([k k+1])=[];
    case 'timeclipmin',
      TimeClipMin=varargin{k+1};
      varargin([k k+1])=[];
    case 'timeclipmax',
      TimeClipMax=varargin{k+1};
      varargin([k k+1])=[];
    case 'asd',
      ASD=varargin{k+1};
      if length(ASD) ~=6
          error('Length of ASD must be 6.')
      end
      ASD=datenum(ASD(1),ASD(2),ASD(3),ASD(4),ASD(5),ASD(6));
      varargin([k k+1])=[];
   otherwise
      k=k+2;
  end;
end;

if length(varargin)<2
   varargin={};
end

disp(['Gregorian Time Offset = ' num2str(ASD)])

[fid,message]=fopen(FortName,'r');
if fid==-1
   error(['Could not open ' FortName ' because ' message])
end

% The header line: RUNDES, RUNID, AGRID 
header=fgets(fid);

% NTRSPE, NSTAE, DT*NSPOOLE, NSPOOLE, IRTYPE 
data=fscanf(fid,'%d %d %f %d %d',[5])';
fgetl(fid);

NTRSPE=data(1);
NSTAE=data(2);
outdt=data(3);
NSPOOLE=data(4);

% determine clipping times
t0=outdt/86400;
tend=outdt*NTRSPE/86400;
tt=t0:outdt/86400:tend;
if isnan(TimeClipMin)
   i0=1;
else
   [a,i0]=min(abs(tt-TimeClipMin));
end
if isnan(TimeClipMax)
   iend=NTRSPE;
else
   [a,iend]=min(abs(tt-TimeClipMax));
end

zeta=NaN*ones(NSTAE,iend-i0+1);
t=NaN*ones(iend-i0+1,1);

jj=0;
ii=0;

if (~Compact)
   while ~feof(fid)
      ii=ii+1;
      temp=fscanf(fid,'%f %d',[1 2]);
      data=fscanf(fid,'%d %f',[2 NSTAE])';
      
      if (ii>=i0) & (ii<=iend) 
         disp(['Storing ' num2str(temp(1)) ])
         jj=jj+1;
         zeta(:,jj)=data(:,2);
         iter(jj)=temp(2);
         t(jj)=temp(1);
      else
         disp(['Skipping ' num2str(temp(1)) ])
      end
      
      fgetl(fid);
   end
else
   while ~feof(fid)
      ii=ii+1;
      temp=fscanf(fid,'%f %d %d %f',[1 4]);
      NWET=temp(3);
      if (ii>=i0) & (ii<=iend) 
         fprintf('Storing iter=%d : Time=%f  : NWET = %d\n',temp(1),temp(1)/86400,NWET)
         jj=jj+1;
         fillval=temp(4);
         t(jj)=temp(1);      
         iter(jj)=temp(2);
         data=fscanf(fid,'%d %f',[2 NWET])';
         zeta(:,jj)=fillval;
         zeta(data(:,1),jj)=data(:,2);
      else
         disp(['Skipping ' num2str(temp(1)) ])
         data=fscanf(fid,'%d %f',[2 NWET])';
      end
      fgetl(fid);
   end
end


if isfinite(ClipVal)
   inan=find(zeta<ClipVal);
   zeta(inan)=NaN;
end

%Fill return structure

D.zeta=zeta;
D.t=t;
D.iter=iter;
D.t=D.t/86400+ASD;
D.nodes=NSTAE;
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
