function D=TidalDatumsFromTimeseries(z,t)
% computes tidal datums and CDFs from s time series, presumed a tidal time
% series

edges=-5:.01:5;  % for PDFs,CDFs

NEdges=length(edges);
CDF=NaN*ones(NEdges,1);
PDF=CDF;

MHHW = NaN*ones(1,1);
MHW  = MHHW; 
MSL  = MHHW; 
MLW  = MHHW; 
MLLW = MHHW; 
RMS  = MHHW;

i=1;

    zt=z;  %  nansum(z);
    
    RMS(i)=rms(zt);
    
    % PDF,CDF
    [n,~] = histc(zt,edges);
    ss=sum(n);
    n=n/ss;
    PDF(:,i)=n;
    CDF(:,i)=cumsum(n)';
  
    MSL(i)=mean(zt);

    dzt=diff(zt);

    inilneg=diff(sign(dzt))==2;
    inilpos=diff(sign(dzt))==-2;

    inilpos=[inilpos(end); inilpos(1:end-1)];
    inilneg=[inilneg(end); inilneg(1:end-1)];

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
D.HAT = max(z);
D.LAT = min(z);

%D.HAT = sum(AdcircHaStruct.AMP(:,TheseFreqs),2);
% D.SpringTide=AdcircHaStruct.AMP(:,idxm2)+AdcircHaStruct.AMP(:,idxs2);
% D.NeapTide=AdcircHaStruct.AMP(:,idxm2)-AdcircHaStruct.AMP(:,idxs2);
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
