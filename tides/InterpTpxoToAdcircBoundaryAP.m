function ReturnThisZeta=InterpTpxoToAdcircBoundaryAP(fgs,varargin)
%function ReturnThisZeta=InterpTpxoToAdcircBoundaryAP(fgs,TpxoModelFileUrl,varargin)
% InterpTpxoToAdcircBoundaryAP(FemGridStruct,TpxoModelFile,...)
% ca=InterpTpxoToAdcircBoundaryAP(fgs,TpxoModelFileUrl,varargin)

TpxoModelFileUrl='http://tds.renci.org:8080/thredds/dodsC/DataLayers/Tides/TPXO/h_tpxo7.2.nc';

%% extract model zeta 
% file=sprintf('%s/%s',TpxoModelPath,TpxoModelFile);
% if ~exist(file)
%    error('TPXO Model file %s not found.',file)
% end

zeta=ncgeodataset(TpxoModelFileUrl);
temp=zeta{'con'};
temp=temp(:,:);

for i=1:length(temp)
   tpxo_constits{i}=deblank(temp(i,:));
end
x=zeta{'lon_z'};x=x(:,:);
y=zeta{'lat_z'};y=y(:,:);
x=x';
y=y';

% report possible constits
%disp(sprintf('Constituents in TPXO model file %s',tpxo_constits))

%% process varargins
ShowPlot=0;
Report=1;
Constits=tpxo_constits;
ReportFileName='tidal_constits.report';

k=1;
while k<length(varargin)
  switch lower(varargin{k})
    case 'showplot'
      ShowPlot=varargin{k+1};
      varargin([k k+1])=[];
    case 'report'
      Report=varargin{k+1};
      varargin([k k+1])=[];
    case 'reportfilename'
      ReportFileName=varargin{k+1};
      varargin([k k+1])=[];
    case 'constits'
      Constits=varargin{k+1};
      if ~iscell(Constits)
         error('Constituent argument must be a cell array.')
      end
      varargin([k k+1])=[];
      for  i=1:length(Constits)
         idx=find(strcmp(Constits{i},tpxo_constits));
         if isempty(idx)
            error('Constituent %s not found in tpxo data file. Terminal.',Constits{i})
         end
      end
    otherwise
      k=k+2;
  end
end

if length(varargin)<2
   varargin={};
end


ReturnThisZeta=[];

%% Accumulate open boundary nodes
xb=NaN*ones(fgs.nopenboundarynodes{1},1);
yb=xb;

for i=1:fgs.nopenboundarynodes{1}
   xb(i)=fgs.x(fgs.ob{1}(i));
   yb(i)=fgs.y(fgs.ob{1}(i));
end
% if the min x (long) in the tpxo grid is positive, assume that 
% the grid is in the range [0 360].  So, put ADCIRC boundary nodes into
% same range.
if min(x(:)>0)
   xb=360 + xb;
end

%% Interpolate constituents
hRe=zeta{'hRe'};
hIm=zeta{'hIm'};
for i=1:length(Constits)
   
   k=find(strcmp(Constits{i},tpxo_constits));

   GlobalZetaAmplitude=squeeze(zeta{'ha'}(k,:,:))';
   
   GlobalZetaPhase=squeeze(zeta{'hp'}(k,:,:))'*pi/180;
   
   iding=find(abs(GlobalZetaAmplitude)==0);
   GlobalZetaAmplitude(iding)=NaN;
   GlobalZetaPhase(iding)=NaN;
      
   GlobalZeta=GlobalZetaAmplitude.*exp(1i*GlobalZetaPhase);
   
   ThisZeta(:,i)=interp2(x,y,GlobalZeta,xb,yb);
    
   % fix NaNs with nearest neighbors in tpxo solution
   inan=find(isnan(ThisZeta(:,i)));
   if ~isempty(inan)
      dxyT=min(diff(unique(x(:))));
      for ii=1:length(inan)
         xx=xb(inan(ii))+5*[-1 1]*dxyT;
         yy=yb(inan(ii))+5*[-1 1]*dxyT;
         idx=find((x>xx(1) & x<xx(2) & y>yy(1) & y<yy(2)) & isfinite(GlobalZeta) );
         d=(x(idx)-xb(inan(ii))).^2 + (y(idx)-yb(inan(ii))).^2;
         [~,idxd]=min(d);
         ThisZeta(inan(ii),i)=GlobalZeta(idx(idxd));
      end   
   end
   
   ThisZetaAmplitude(:,i)=abs(ThisZeta(:,i));
   ThisZetaPhase(:,i)=mod(angle(ThisZeta(:,i))*180/pi,360);
   
   if ShowPlot
      if (i==1)
         fgs.x=fgs.x+360;
         dxr=(max(xb)-min(xb))/10;
         dyr=(max(yb)-min(yb))/10;
         x0=min(xb);
         x1=max(xb);
         y0=min(yb);
         y1=max(yb);

      end

      ca=[min(ThisZetaAmplitude(:,i)) max(ThisZetaAmplitude(:,i))];
      cp=[min(ThisZetaPhase(:,i))     max(ThisZetaPhase(:,i))];

      figure
        %subplot(121)
        pcolor(x,y,GlobalZetaAmplitude)
	     shading flat
	     axis('equal')
	     set(gca,'CLim',ca);
	     colorbar
	     axis([x0-dxr x1+dxr y0-dyr y1+dyr])
	     line(xb,yb,'Marker','*','Color','k')
	     plotbnd(fgs,'Color','k')
	     numscal(xb(1:1:end),yb(1:1:end),ThisZetaAmplitude(1:1:end,i),'FontSize',18)
	     title({upper(Constits{i}),'Amplitude'},'Fontsize',13)
       
        figure
        %subplot(122)
	     pcolor(x,y,GlobalZetaPhase*180/pi)
	     shading flat
	     axis('equal')
	     set(gca,'CLim',cp);
	     colorbar
	     axis([x0-dxr x1+dxr y0-dyr y1+dyr])
	     line(xb,yb,'Marker','*','Color','k')
	     plotbnd(fgs,'Color','k')
	     numscal(xb(1:10:end),yb(1:10:end),ThisZetaPhase(1:10:end,i),'FontSize',18)
	     title({upper(Constits{i}),'Phase'},'Fontsize',13) 
   end  
   

   
   
end

%% 
if Report
   fid=fopen(ReportFileName,'w');
   for i=1:length(Constits)
      fprintf(fid,'%s\n',upper(Constits{i}));
      out=[ThisZetaAmplitude(:,i),ThisZetaPhase(:,i)];
      fprintf(fid,'%f %f\n',out');
   end
end

ReturnThisZeta.amp=ThisZetaAmplitude;
ReturnThisZeta.pha=ThisZetaPhase;
ReturnThisZeta.Constits=Constits;   

return

