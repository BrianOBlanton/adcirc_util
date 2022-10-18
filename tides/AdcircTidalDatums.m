function D=AdcircTidalDatums(AdcircHaStruct,LenYears)
% computes tidal datums and CDFs from ADICRC harmonic analysis struct
% D=AdcircTidalDatums(AdcircHaStruct)

% Constits={'M2','S2','O1','K1'};
% %[TheseConstits.c,TheseConstits.ia,TheseConstits.ib]=intersect(AdcircHaStruct.PERNAMES,Constits);
% 
% for i=1:length(Constits)
%     idx(i)=find(strcmp(AdcircHaStruct.PERNAMES,Constits(i)));
% end

%BlockSize=10000;

if ~exist('LenYears','var'),LenYears=1;end

NStations=size(AdcircHaStruct.AMP,1);
%NFreqs=length(AdcircHaStruct.FREQ);
% freq in cph, period in hrs
Freqs=AdcircHaStruct.FREQ*3600; 
%Periods=2*pi./Freqs;D=ans;
edges=-5:.01:5;  % for PDFs,CDFs

idxz0=find(strcmp(AdcircHaStruct.PERNAMES,'STEADY'));
AdcircHaStruct.PERNAMES(idxz0)=[];
AdcircHaStruct.PHA(:,idxz0)=[];
AdcircHaStruct.AMP(:,idxz0)=[];
AdcircHaStruct.FREQ(idxz0)=[];


%           'M2' 'S2' 'N2' 'K2' 'O1' 'K1' 'P1' 'Q1'
UseTheseConstits={'M2' 'S2' 'N2' 'K2' 'O1' 'K1' 'P1' 'Q1'};
%UseTheseConstits={'M2' 'S2' 'N2' 'O1' 'K1' 'Q1'};
[~,TheseFreqs]=ismember(UseTheseConstits,AdcircHaStruct.PERNAMES);
%TheseFreqs=[2:23];
%TheseFreqs

NEdges=length(edges);
CDF=NaN*ones(NStations,NEdges);
PDF=CDF;

idxm2=find(strcmp(AdcircHaStruct.PERNAMES,'M2'));
idxs2=find(strcmp(AdcircHaStruct.PERNAMES,'S2'));

MHHW = NaN*ones(NStations,1);
MHW  = MHHW; 
MSL  = MHHW; 
MLW  = MHHW; 
MLLW = MHHW; 
RMS  = MHHW;

inan=AdcircHaStruct.AMP(:,idxm2)<1e-2;  % locations at which m2 amp < tol
AdcircHaStruct.AMP(inan,:)=NaN; % set these locations to NaN
AdcircHaStruct.PHA(inan,:)=NaN;

inan=isnan(AdcircHaStruct.PHA);
AdcircHaStruct.PHA(inan)=0;

sum_amp=nansum(AdcircHaStruct.AMP,2);
inan=sum_amp<1e-2;
AdcircHaStruct.AMP(inan,:)=NaN; % set these locations to NaN
AdcircHaStruct.PHA(inan,:)=NaN;


% ti=0;            % starting tidal cycle
% te=706;          % ending tidal cycle ; 706 is 365.3550 days; 13401 is 6.9350e+03 days (19 years)
% dt=12.42/50;     % hours
% ti=ti*Periods(idxm2);
% te=te*Periods(idxm2);

ti=0;               % starting time
te=LenYears*365.25; % ending time in days
dt=.25;             % hours
ti=ti*24;
te=te*24;

t=(ti:dt:te)';
t(end)=[];

amps=AdcircHaStruct.AMP(:,TheseFreqs); 
phas=AdcircHaStruct.PHA(:,TheseFreqs)*pi/180; 
freqs=Freqs(TheseFreqs);
M2amp=AdcircHaStruct.AMP(:,idxm2); 
M2pha=AdcircHaStruct.PHA(:,idxm2); 

BlockSize=1000;

for i=1:NStations
    
    if mod(i,BlockSize)==0
        fprintf('Pct Complete = %5.1f\n',100*i/NStations)
    end
     
    if isnan(M2amp(i)) || isnan(M2pha(i))
        continue
    end
    
    z=NaN*ones(length(TheseFreqs),length(t));
    
    for j=1:length(TheseFreqs)
        fj=freqs(j);
        z(j,:)=amps(i,j)*cos(fj*t+phas(i,j)*pi/180);
    end
    zt=nansum(z);
    
    RMS(i)=rms(zt);
    
    % PDF,CDF
    [n,~] = histc(zt,edges);
    ss=sum(n);
    n=n/ss;
    PDF(i,:)=n;
    CDF(i,:)=cumsum(n);
  
    MSL(i)=mean(zt);

    dzt=diff(zt);

    inilneg=diff(sign(dzt))==2;
    inilpos=diff(sign(dzt))==-2;

    inilpos=[inilpos(end) inilpos(1:end-1)];
    inilneg=[inilneg(end) inilneg(1:end-1)];

    highs=zt(inilpos);
    lows=zt(inilneg);

    % mean high water
    MHW(i)=mean(highs);
 
    % mean low water
    MLW(i)=mean(lows);

    if rem(length(highs),2)==1
       highs(end)=[];
    end
    if rem(length(lows),2)==1
       lows(end)=[];
    end

    if highs(1)>highs(2)
       temp=reshape(highs,2,length(highs)/2);
    else
       highs([1 end])=highs([end 1]);
       temp=reshape(highs,2,length(highs)/2);
    end
    temp2=sort(temp,1,'descend');
    MHHW(i)=mean(temp2(1,:));

    if lows(1)<lows(2)
       temp=reshape(lows,2,length(lows)/2);
    else
       lows([1 end])=lows([end 1]);
       temp=reshape(lows,2,length(lows)/2);
    end
    temp2=sort(temp,1,'ascend');
    MLLW(i)=mean(temp2(1,:));

end

D.PDF=PDF;
D.CDF=CDF;
D.edges=edges;
D.MHHW=MHHW;
D.MHW=MHW;
D.MSL=MSL;
D.MLW=MLW;
D.MLLW=MLLW;

D.MTL = mean([D.MHW D.MLW],2);
D.DTL = mean([D.MHHW D.MLLW],2);
D.MN  = D.MHW-D.MLW;
D.GT  = D.MHHW-D.MLLW;
D.HAT = sum(AdcircHaStruct.AMP(:,TheseFreqs),2);
D.SpringTide=AdcircHaStruct.AMP(:,idxm2)+AdcircHaStruct.AMP(:,idxs2);
D.NeapTide=AdcircHaStruct.AMP(:,idxm2)-AdcircHaStruct.AMP(:,idxs2);
D.RMS=RMS;


% D.MSL(abs(D.MSL)<.0000001) = NaN;
% D.MLLW(D.MLLW>-.0001)      = NaN;
% D.MHHW(D.MHHW<.0001)       = NaN;
% D.MLW(D.MLW>-.0001)        = NaN;
% D.MHW(D.MHW<.0001)         = NaN;

% D.HAT(D.LAT<.0001)         = NaN;

return


clf
plot(t/24,zt,'-')
line(t(inilneg)/24,lows,'Marker','x','Color','g','LineStyle','none')
line(t(inilpos)/24,highs,'Marker','x','Color','r','LineStyle','none')
   %hline(mhw,'Color','r');
       %hline(mlw,'Color','g');
    %hline(mhhw,'Color','r','Linewidth',2);
        %hline(mllw,'Color','g','Linewidth',2);
