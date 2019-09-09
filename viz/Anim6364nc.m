function Anim6364nc(varargin)
% Anim6364nc generate an animatiom of an ADCIRC 63,64 netCDF file pair
% Defaults: 
%     nc63 - will open fort.63.nc as an ncgeodataset
%     nc64 - will open fort.64.nc as an ncgeodataset
%     fgs  - will be extracted from nc63
%
% P/V pairs:
%     Grid
%     nc63
%     nc64
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
%     ImageWriteDir
%     FrameBaseName - base of image output file name (def='frame')
%     ManualAdvance
%     Url 

%     StartingTime - datenum of starttime

% if nargin==0
%    disp('setting defaults:')
%    nc=ncgeodataset('fort.63.nc');
%    nc64=ncgeodataset('fort.64.nc');
%    g=ExtractGrid(nc);
% end

p.FrameBaseName='frame';
p.ScriptToAdd='none';
p.ImageResolution='-r100';
p.ImageWriteDir='./';
p.AxisLims=[];  
p.Title={};     
p.IterStart=1; 
p.IterStride=1;
p.IterStop=-1;  
p.ColorMin=NaN;  
p.ColorMax=NaN;  
p.ColorMap=jet(32);  
%StartingTime=0;
p.VectorStride=1;
p.VectorScaleFac=1;
p.VectorColor='k';
p.ScaleXor=[];
p.ScaleYor=[];
p.ScaleLabel='no scale';
p.ManualAdvance=false;
p.Grid='';
p.nc63='fort.63.nc';
p.nc64='fort.64.nc';
p.Url='';

p=parse_pv_pairs(p,varargin);

if isempty(p.Url)
    if exist('p.nc63','file')
        p.nc63=ncgeodataset(p.nc63);
        p.nc64=ncgeodataset(p.nc64);
    else
        error('fort.63,64.nc does not exist locally. Either cd to a directory, or use Url.')
    end
else
    try
        p.nc63=ncgeodataset([p.Url '/fort.63.nc']);
        p.nc64=ncgeodataset([p.Url '/fort.64.nc']);
    catch ME
        disp(['ID: ' ME.identifier])
        throw(ME)
    end
end

if isempty(p.Grid)
    p.Grid=ExtractGrid(p.nc63);
end

time=p.nc63.time('time');
time=datetime(datevec(time));
nTimes=size(time,1); 
t=time;

if p.IterStop==-1,p.IterStop=nTimes;end
if p.IterStart==-1,p.IterStart=nTimes;p.IterStop=nTimes;end

if (p.IterStart<1 || p.IterStart>nTimes)
   error('IterStart must be between %d and %d',1,nTimes)
elseif (p.IterStop<1 || p.IterStop>nTimes)
   error('IterStop must be between %d and %d',1,nTimes)
elseif (p.IterStop<p.IterStart)
   error('IterStop must be greater than or equal to IterStart')
end

if p.IterStride<1
    error('Stride must be greater than 1.')
end

fprintf('Starting iteration =  %d\n',p.IterStart)
fprintf('Iteration Stride   =  %d\n',p.IterStride)
fprintf('Stopping iteration =  %d\n',p.IterStop)

h=[];
hu0=[];
%hz0=[];
%hst=[];
%axx=[];

if ~iscell(p.Title)
   error('Title to Anim63 must be a cell array')
end

tlen=length(p.Title);

% set up figure
if (1)
    figure
    %drawelems(g,'Color',[1 1 1]*.7,'Linewidth',.25);
    plotbnd(p.Grid,'LineWidth',.2);
    hc=lcontour(p.Grid,'z',0,'Color','k','LineWidth',.2);
    if isgraphics(hc)
        set_height(hc,1);
    end
    plotcoast('states')
    grid on
    grid minor
    %hc=lcontour(g,'z',[2:8],'Color','r','LineWidth',.2);
    %set_height(hc,1);
    %hc=lcontour(g,'z',-[2:10],'Color','b','LineWidth',.2);
    %set_height(hc,1);
end
axis('equal')

if isempty(p.AxisLims)
    minX=min(p.Grid.x);maxX=max(p.Grid.x);
    minY=min(p.Grid.y);maxY=max(p.Grid.y);
else
    minX=p.AxisLims(1);
    maxX=p.AxisLims(2);
    minY=p.AxisLims(3);
    maxY=p.AxisLims(4);
end
axis([minX maxX minY maxY])

% get grid indices in box;
InViewingBox=find(p.Grid.x>minX & p.Grid.x<maxX & ...
                  p.Grid.y>minY & p.Grid.y<maxY);

if exist(p.ScriptToAdd,'file')
   eval(p.ScriptToAdd)
end

if ~(isnan(p.ColorMin) && isnan(p.ColorMax))    
    set(gca,'CLim',[p.ColorMin p.ColorMax])
end
colormap(p.ColorMap)
colorbar

% generate frames
for i=p.IterStart:p.IterStride:p.IterStop
   
   if isnat(t(i)),break,end
   %zz=D.zeta(:,i);
   uu=p.nc64{'u-vel'}(i,:)';
   vv=p.nc64{'v-vel'}(i,:)';
   zz=p.nc63{'zeta'}(i,:)';
   zz(zz<-1000)=NaN;
   if ~isempty(h),delete(h),end   
   if ~isempty(hu0),delete(hu0),end
   %if ~isempty(hz0),delete(hz0),end
   %if ~isempty(hst),delete(hst);delete(axx),end
   h=colormesh2d(p.Grid,zz);
   %hz0=lcontour(g,zz,0,'Color','b');
   hu0=vecplot(p.Grid.x,p.Grid.y,uu,vv,...
       'ScaleLabel',p.ScaleLabel,...
       'ScaleXor',p.ScaleXor,...
       'ScaleYor',p.ScaleYor,...
       'ScaleFac',p.VectorScaleFac,...
       'Stride',p.VectorStride,...
       'Color',p.VectorColor);
   
   if ~(isnan(p.ColorMin) && isnan(p.ColorMax))    
       set(gca,'CLim',[p.ColorMin p.ColorMax])
   else
       cmin=min(zz(InViewingBox));
       cmax=max(zz(InViewingBox));
       set(gca,'CLim',[cmin cmax])
   end
   %set(gca,'CLim',[-1 1]*max(zz))
   %shading flat
   %set_height(h,1)
   p.Title{tlen+1}=sprintf('%s  //  %04d',datestr(t(i),0),i);
   title(p.Title,'FontSize',14);
   %[hst,axx]=stamp_right(datestr(t(i)));
   drawnow

   fnamebase=sprintf('%s_%03d',p.FrameBaseName,i);
   fprintf('Printing %s\n',fnamebase)
   export_fig('-png',p.ImageResolution,sprintf('%s/%s.png',p.ImageWriteDir,fnamebase));
%   print('-dpng',ImageResolution,sprintf('%s.png',fnamebase));
   %eval(sprintf('!/usr/local/bin/convert %s.png %s.gif',fnamebase,fnamebase));
    if p.ManualAdvance
        pause
    end
end

