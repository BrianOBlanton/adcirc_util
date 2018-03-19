%function TidalSolutionComparision(f15,f51,MainStationList,ThisStationList)
function TidalSolutionComparison(f14,ThisStationList,ThisConstitList)
% Call as: TidalSolutionComparison(f14,ThisStationList,ThisConstitList);

NosInit;

% Defaults
%f14='fort.14';
f15='fort.15';
f53='fort.51';
MainStationList='station.list';
Constits={'M2' 'N2' 'S2' 'K2' 'O1' 'K1' 'P1' 'Q1'};
PlotConstits={'M2' 'K1'};

if nargin==0
    error('Need station list cell array as input.')
end

if ~iscell(ThisStationList)
    error('ThisStationList must be a cell array')
end

% constite to process; must be a subset of constits analyzed by ADCIRC
if ~exist('ThisConstitList')
else
    Constits=ThisConstitList;
end

%f14=grd_to_opnml;
f15=read_adcirc_fort15('fort.15',f14.nopennodes{1});
f51=read_adcirc_fort53('FileName',f53);

% get ADCIRC station list
[AdcircStationList.x,AdcircStationList.y,AdcircStationList.id]=...
	textread('../../station.list','%f%f%s%*[^\n]');

%%
% get NOS tides at these stations
NosTides=LoadNosTides(ThisStationList);

% find stations in Master Station List
[TheseStations.c,TheseStations.ia,TheseStations.ib]=intersect(AdcircStationList.id,ThisStationList);
if length(ThisStationList) ~= length(TheseStations.c)
    disp('Some stations not found in master station list...')
    disp(setdiff(ThisStationList,TheseStations.c))
end
[TheseConstits.InAdc.c,TheseConstits.InAdc.ia,TheseConstits.InAdc.ib]=intersect(f51.PERNAMES,Constits);
[TheseConstits.InNos.c,TheseConstits.InNos.ia,TheseConstits.InNos.ib]=intersect(NosTides.Constituent(:,1),Constits);

% keep only tides in Constits cell array
NosTides.Constituent = NosTides.Constituent(TheseConstits.InNos.ia,:);
NosTides.Amplitude   = NosTides.Amplitude(TheseConstits.InNos.ia,:);
NosTides.Phase       = NosTides.Phase(TheseConstits.InNos.ia,:);
NosTides.Speed       = NosTides.Speed(TheseConstits.InNos.ia,:);
NosTides.freq        = NosTides.freq(TheseConstits.InNos.ia,:);
NosTides.ampdata     = NosTides.ampdata(TheseConstits.InNos.ia,:);
NosTides.phadata     = NosTides.phadata(TheseConstits.InNos.ia,:);

AdcTides.Constituent = f51.PERNAMES(TheseConstits.InAdc.ia)';
AdcTides.freq        = f51.FREQ(TheseConstits.InAdc.ia)';
AdcTides.Amplitude   = f51.AMP(TheseStations.ia,TheseConstits.InAdc.ia)';
AdcTides.Phase       = f51.PHA(TheseStations.ia,TheseConstits.InAdc.ia)';

CmplxAmp.Adc=AdcTides.Amplitude.*exp(sqrt(-1)*AdcTides.Phase*pi/180);
CmplxAmp.Nos=NosTides.Amplitude.*exp(sqrt(-1)*NosTides.Phase*pi/180);
CmplxAmp.Dif=CmplxAmp.Nos-CmplxAmp.Adc;
CmplxError=abs(CmplxAmp.Dif);

%% default plots
for i=1:length(PlotConstits)
    idxnos=strmatch(PlotConstits{i},NosTides.Constituent(:,1));
    idxadcirc=strmatch(PlotConstits{i},AdcTides.Constituent(:,1));
    figure
    o=NosTides.Amplitude(idxnos,:);
    a=AdcTides.Amplitude(idxadcirc,:);
   subplot(121)
      plot(o,a,'Marker','*','LineStyle','none')
      line([0 1],[0 1])
      axeq
      grid
      title(sprintf('%s Amp [m]',AdcTides.Constituent{idxadcirc}))
   subplot(122)
      o=NosTides.Phase(idxnos,:);
      a=AdcTides.Phase(idxadcirc,:);
      %if idx==3
      %    o(o>180)=o(o>180)-360;
      %    a(a>180)=a(a>180)-360;
      %end
      plot(o,a,'Marker','*','LineStyle','none')
      %if idx==3
      %    line([-10 20],[-10 20])
      %    line([-5 20],[-10 15])
      %    line([-10 15],[-5 20])
      %    axis([-10 20 -10 20])
      %else
      %    line([0 360],[0 360])
      %end
      axeq
      grid
      title(sprintf('%s Pha [deg]',AdcTides.Constituent{idxadcirc}))
end

%%
[m,n]=size(NosTides.Amplitude);
da=100*(1-NosTides.Amplitude./AdcTides.Amplitude);
dp=NosTides.Phase-AdcTides.Phase;
fmtstr  = ['Const  ' repmat('%-7s', [1 n]) '\n'];
fmtstr2 = ['%-4s'    repmat('%7.2f',[1 n]) '\n'];
fmtstr2p= ['%-4s'    repmat('%7.1f',[1 n]) '\n'];

fid=fopen('aamp.txt','w');
fprintf(fid,fmtstr,NosTides.StationNameAbbr{:});
for i=1:length(NosTides.Amplitude)
   fprintf(fid,fmtstr2,AdcTides.Constituent{i},AdcTides.Amplitude(i,:));
end
fid=fopen('namp.txt','w');
fprintf(fid,fmtstr,NosTides.StationNameAbbr{:});
for i=1:length(NosTides.Amplitude)
   fprintf(fid,fmtstr2,NosTides.Constituent{i},NosTides.Amplitude(i,:));
end
fid=fopen('damp.txt','w');
fprintf(fid,fmtstr,NosTides.StationNameAbbr{:});
for i=1:length(NosTides.Amplitude)
   fprintf(fid,fmtstr2,NosTides.Constituent{i},da(i,:));
end

fid=fopen('apha.txt','w');
fprintf(fid,fmtstr,NosTides.StationNameAbbr{:});
for i=1:length(NosTides.Amplitude)
   fprintf(fid,fmtstr2p,AdcTides.Constituent{i},AdcTides.Phase(i,:));
end
fid=fopen('npha.txt','w');
fprintf(fid,fmtstr,NosTides.StationNameAbbr{:});
for i=1:length(NosTides.Amplitude)
   fprintf(fid,fmtstr2p,NosTides.Constituent{i},NosTides.Phase(i,:));
end
fid=fopen('dpha.txt','w');
fprintf(fid,fmtstr,NosTides.StationNameAbbr{:});
for i=1:length(NosTides.Amplitude)
   fprintf(fid,fmtstr2p,NosTides.Constituent{i},dp(i,:));
end
fclose(fid);

