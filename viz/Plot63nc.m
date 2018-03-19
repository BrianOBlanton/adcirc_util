function AxesHandle=Plot63nc(g,nc,varargin)
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




if nargin==0
   disp('Plot63nc(g,nc) OR:')
   disp('Plot63nc(fgs,nc,p1,v1,p2,v2,...)')
   return
end

% set defaults
% FrameBaseName='frame';
ScriptToAdd='none';
% ImageResolution='-r200';
AxisLims=[];  
Title={''};     
ColorMin=NaN;  
ColorMax=NaN;  
ColorMap=jet(32);  
StartingTime=0;
Iter=NaN;
SetUpAxes=true;
AxesHandle=NaN;

% Strip off propertyname/value pairs in varargin not related to
% "line" object properties.
k=1;
while k<length(varargin),
  switch lower(varargin{k}),
    case 'startingtime',
      StartingTime=varargin{k+1};
      varargin([k k+1])=[];
    case 'framebasename',
      FrameBaseName=varargin{k+1};
      varargin([k k+1])=[];
    case 'scripttoadd',
      ScriptToAdd=varargin{k+1};
      varargin([k k+1])=[];
    case 'axislims',
      AxisLims=varargin{k+1};
      varargin([k k+1])=[];
    case 'iter',
      Iter=varargin{k+1};
      varargin([k k+1])=[];
    case 'title',
      Title=varargin{k+1};
      varargin([k k+1])=[];
    case 'colormin',
      ColorMin=varargin{k+1};
      varargin([k k+1])=[];
    case 'colormax',
      ColorMax=varargin{k+1};
      varargin([k k+1])=[];
    case 'colormap',
      ColorMap=varargin{k+1};
      varargin([k k+1])=[];
    case 'setupaxes',
      SetUpAxes=varargin{k+1};
      varargin([k k+1])=[];
    case 'axeshandle',
      AxesHandle=varargin{k+1};
      varargin([k k+1])=[];
    otherwise
      k=k+2;
  end;
end;

if length(varargin)<2
   varargin={};
end

nTimes=ncsize(nc{'time'});
t=StartingTime+nc{'time'}(:)/24;

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


if ~iscell(Title)
   error('Title to Plot63nc must be a cell array')
end

tlen=length(Title);

if exist(ScriptToAdd,'file')
   eval(ScriptToAdd)
end

if ~(isnan(ColorMin) && isnan(ColorMax))    
    set(gca,'CLim',[ColorMin ColorMax])
end
colormap(ColorMap)
colorbar

% generate frame
   
   ff=get(gcf,'Position');
   axes(AxesHandle)
   %set(gcf,'Position',ff)
   
   if isnan(t(Iter)),return,end
   %zz=D.zeta(:,i);
   zz=nc{'zeta'}(Iter,:)';
   zz(zz<-1000)=NaN;
   %if ~isempty(h),delete(h),end   
   %if ~isempty(hz0),delete(hz0),end
   %if ~isempty(hst),delete(hst);delete(axx),end
   delete(findobj(AxesHandle,'Tag','colorsurf'));
   h=colormesh2d(g,zz);
   %hz0=lcontour(g,zz,0,'Color','b');
   if ~(isnan(ColorMin) & isnan(ColorMax))    
       set(gca,'CLim',[ColorMin ColorMax])
   else
       cmin=min(zz(InViewingBox));
       cmax=max(zz(InViewingBox));
       set(gca,'CLim',[cmin cmax])
   end
   %set(gca,'CLim',[-1 1]*max(zz))
   %shading flat
   %set_height(h,1)
   Title{tlen+1}=datestr(t(Iter),0);
   title(Title,'FontSize',14);
   %[hst,axx]=stamp_right(datestr(t(i)));
   drawnow


%end
