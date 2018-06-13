function Anim63nc(g,nc,varargin)
% Anim63nc generate an animatiom of an ADCIRC 63 file in netCDF format
% Required inputs: 
%     fgs - grid structure on which the 63 file was computed 
%     nc  - netcdf object from "netcdf" pointing to a fort.63.nc file
% 
% P/V pairs:
%     AxisLims   - plot axis limits, as in the axis command (def=grid lims)
%     Title      - title string as a cell array (def={''})
%     IterStart  - starting iteration number (def=1)
%     IterStride - iteration stride (skip) (def=1)
%     IterStop   - stopping iteration number (def=end)
%     ColorMin   - min pressure to clip below (def=950)
%     ColorMax   - max pressure to clip above (def=1030)
%     ColorMap   - colormap to use (def=jet(32))
%     ScriptToAdd - script that defines plot overlays (def='none')
%     ImageResolution - (def='-r200';)
%     FrameBaseName - base of image output file name (def='frame')
%     StartingTime - datenum of starttime

if nargin==0
   disp('Anim63nc(g,nc) OR:')
   disp('Anim63nc(fgs,nc,p1,v1,p2,v2,...)')
   return
end

FrameBaseName='frame';
ScriptToAdd='none';
ImageResolution='-r200';
AxisLims=[];  
Title={''};     
IterStart=1; 
IterStride=1;
IterStop=-1;  
ColorMin=NaN;  
ColorMax=NaN;  
ColorMap=jet(32);  
StartingTime=0;

% Strip off propertyname/value pairs in varargin not related to
% "line" object properties.
k=1;
while k<length(varargin)
  switch lower(varargin{k})
    case 'startingtime'
      StartingTime=varargin{k+1};
      varargin([k k+1])=[];
    case 'framebasename'
      FrameBaseName=varargin{k+1};
      varargin([k k+1])=[];
    case 'scripttoadd'
      ScriptToAdd=varargin{k+1};
      varargin([k k+1])=[];
    case 'axislims'
      AxisLims=varargin{k+1};
      varargin([k k+1])=[];
    case 'title'
      Title=varargin{k+1};
      varargin([k k+1])=[];
    case 'iterstart'
      IterStart=varargin{k+1};
      varargin([k k+1])=[];
    case 'iterstride'
      IterStride=varargin{k+1};
      varargin([k k+1])=[];
    case 'iterstop'
      IterStop=varargin{k+1};
      varargin([k k+1])=[];
    case 'colormin'
      ColorMin=varargin{k+1};
      varargin([k k+1])=[];
    case 'colormax'
      ColorMax=varargin{k+1};
      varargin([k k+1])=[];
    case 'colormap'
      ColorMap=varargin{k+1};
      varargin([k k+1])=[];
    otherwise
      k=k+2;
  end;
end;

if length(varargin)<2
   varargin={};
end

nTimes=nc.size('time');
t=StartingTime+nc{'time'}(:)/24;
%t=D.t/86400;

if IterStop==-1,IterStop=nTimes;end

if (IterStart<1 || IterStart>nTimes)
   error('IterStart must be between %d and %d',1,nTimes)
   IterStart=1;
elseif (IterStop<1 || IterStop>nTimes)
   error('IterStop must be between %d and %d',1,nTimes)
elseif (IterStop<IterStart)
       error('IterStop must be greater than or equal to IterStart')
end

if IterStride<1
    error('Stride must be greater than 1.')
end

fprintf('Starting iteration =  %d\n',IterStart)
fprintf('Stride =  %d\n',IterStride)
fprintf('Stopping iteration =  %d\n',IterStop)

h=[];
hz0=[];
hst=[];
axx=[];

if ~iscell(Title)
   error('Title to Anim63 must be a cell array')
end

tlen=length(Title);

% set up figure
if (1)
   figure
   %drawelems(g,'Color',[1 1 1]*.7,'Linewidth',.25);
   plotbnd(g,'LineWidth',.2);
   hc=lcontour(g,'z',0,'Color','k','LineWidth',.2);
   set_height(hc,1);
   plotcoast('states')
   grid
   %hc=lcontour(g,'z',[2:8],'Color','r','LineWidth',.2);
   %set_height(hc,1);
   %hc=lcontour(g,'z',-[2:10],'Color','b','LineWidth',.2);
   %set_height(hc,1);
end
axis('equal')

if isempty(AxisLims)
    minX=min(g.x);maxX=max(g.x);
    minY=min(g.y);maxY=max(g.y);
else
    minX=AxisLims(1);
    maxX=AxisLims(2);
    minY=AxisLims(3);
    maxY=AxisLims(4);
end
axis([minX maxX minY maxY])
% get grid indices in box;
InViewingBox=find(g.x>minX & g.x<maxX & g.y>minY & g.y<maxY);

if exist(ScriptToAdd,'file')
   eval(ScriptToAdd)
end

if ~(isnan(ColorMin) && isnan(ColorMax))    
    set(gca,'CLim',[ColorMin ColorMax])
end
colormap(ColorMap)
colorbar

% generate frames
for i=IterStart:IterStride:IterStop
   
   if isnan(t(i)),break,end
   %zz=D.zeta(:,i);
   zz=nc{'zeta'}(i,:)';
   zz(zz<-1000)=NaN;
   if ~isempty(h),delete(h),end   
   if ~isempty(hz0),delete(hz0),end
   %if ~isempty(hst),delete(hst);delete(axx),end
   h=colormesh2d(g,zz);
   %hz0=lcontour(g,zz,0,'Color','b');
   if ~(isnan(ColorMin) && isnan(ColorMax))    
       set(gca,'CLim',[ColorMin ColorMax])
   else
       cmin=min(zz(InViewingBox));
       cmax=max(zz(InViewingBox));
       set(gca,'CLim',[cmin cmax])
   end
   %set(gca,'CLim',[-1 1]*max(zz))
   %shading flat
   %set_height(h,1)
   Title{tlen+1}=datestr(t(i),0);
   title(Title,'FontSize',14);
   %[hst,axx]=stamp_right(datestr(t(i)));
   drawnow

   fnamebase=sprintf('%s_%03d',FrameBaseName,i);
   fprintf('Printing %s\n',fnamebase)
   print('-dpng',ImageResolution,sprintf('%s.png',fnamebase));
   %eval(sprintf('!/usr/local/bin/convert %s.png %s.gif',fnamebase,fnamebase));
end
