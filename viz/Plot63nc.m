function h=Plot63nc(g,nc,varargin)
% Plot63nc draw a timestep of an ADCIRC 63 file in netCDF format
% Required inputs: 
%     fgs - grid structure on which the 63 file was computed 
%     nc  - netcdf object from "netcdf" pointing to a fort.63.nc file
% 
% P/V pairs:
%     AxisLims   - plot axis limits, as in the axis command (def=grid lims)
%     Title      - title string as a cell array (def={''})
%     ColorMin   - min pressure to clip below (def=950)
%     ColorMax   - max pressure to clip above (def=1030)
%     ColorMap   - colormap to use (def=jet(32))
%     ScriptToAdd - script that defines plot overlays (def='none')
%     ImageResolution - (def='-r200';)
%     FrameBaseName - base of image output file name (def='frame')
%     StartingTime - datenum of starttime

% set defaults
% p.FrameBaseName='frame';
p.ScriptToAdd='none';
% p.ImageResolution='-r200';
p.AxisLims=[];  
p.Title={''};     
p.ColorMin=NaN;  
p.ColorMax=NaN;  
p.ColorMap=jet(32);  
p.StartingTime=0;
p.Iter=-1;
p.SetUpAxes=true;
p.AxesHandle=gca;
p.filename='fort.63.nc';
p.Grid='';

p=parse_pv_pairs(p,varargin);

if nargin==0
%   disp('Plot63nc(g,nc) OR:')
%   disp('Plot63nc(fgs,nc,p1,v1,p2,v2,...)')
%   return
    if ~exist(p.filename,'file')
        error('If no arguments, fort.63.nc must locally exist.')
    end   
end

try 
    nc=ncgeodataset(p.filename);
catch
    error('nc file %s not found.',p.filename);
end

if isempty(p.Grid)
    g=ExtractGrid(nc);
else
    g=p.Grid;
end

time=nc.time('time');
time=datetime(datevec(time));

nTimes=size(time,1); 
%t=p.StartingTime+time;
t=time;

% if (IterStart<1 || IterStart>nTimes)
%    error('IterStart must be between %d and %d',1,nTimes)
%    IterStart=1;
% elseif (IterStop<1 || IterStop>nTimes)
%    error('IterStop must be between %d and %d',1,nTimes)
% elseif (IterStop<IterStart)
%        error('IterStop must be greater than or equal to IterStart')
% end

% if IterStride<1
%     error('Stride must be greater than 1.')
% end

% disp(sprintf('Starting iteration =  %d',IterStart))
% disp(sprintf('Stride =  %d',IterStride))
% disp(sprintf('Stopping iteration =  %d',IterStop))

h=[];
hz0=[];
hst=[];
axx=[];

% if ~iscell(p.Title)
%     error('Title to Plot63nc must be a cell array')
% end
% tlen=length(p.Title);

if exist(p.ScriptToAdd,'file')
   eval(p.ScriptToAdd)
end

if ~(isnan(p.ColorMin) && isnan(p.ColorMax))    
    set(gca,'CLim',[p.ColorMin p.ColorMax])
end
colormap(p.ColorMap)
colorbar

if p.Iter==-1
    p.Iter=nTimes;
end

% generate frame
   
   ff=get(gcf,'Position');
   %axes(p.AxesHandle)
   %set(gcf,'Position',ff)
   
   if isnat(t(p.Iter)),return,end
   %zz=D.zeta(:,i);
   zz=nc{'zeta'}(p.Iter,:)';
   zz(zz<-1000)=NaN;
   %if ~isempty(h),delete(h),end   
   %if ~isempty(hz0),delete(hz0),end
   %if ~isempty(hst),delete(hst);delete(axx),end
   delete(findobj(p.AxesHandle,'Tag','colorsurf'));
   h=colormesh2d(g,zz);
   %hz0=lcontour(g,zz,0,'Color','b');
   if ~(isnan(p.ColorMin) && isnan(p.ColorMax))    
       set(gca,'CLim',[p.ColorMin p.ColorMax])
%   else
%       cmin=min(zz(InViewingBox));
%       cmax=max(zz(InViewingBox));
%       set(gca,'CLim',[cmin cmax])
   end
   %set(gca,'CLim',[-1 1]*max(zz))
   %shading flat
   %set_height(h,1)
   %Title{tlen+1}=datestr(t(Iter),0);
   title(p.Title,'FontSize',14);
   %[hst,axx]=stamp_right(datestr(t(i)));
   drawnow


%end
