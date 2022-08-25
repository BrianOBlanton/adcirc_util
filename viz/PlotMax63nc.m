function [h,g]=PlotMax63nc(varargin)
% [h,g]=PlotMax63nc(varargin)
%
% Plot63nc draw an ADCIRC maxele 63 netCDF file
% Required inputs: None.  The default behavior expects a maxele.63.nc file
%                  in the local directory.        
%
% P/V pairs:
%     AxisLims        - plot axis limits, as in the axis command (def=grid lims)
%     ColorMin        - min pressure to clip below (def=950)
%     ColorMax        - max pressure to clip above (def=1030)
%     ColorMap        - colormap to use (def=jet(32))
%     Grid            - ADCIRC grid to use, instead of extracting from nc object
%     ImageResolution - (def='-r200';)
%     FrameBaseName   - base of image output file name (def='frame')
%     ScriptToAdd     - script that defines plot overlays (def='none')
%     StartingTime    - datenum of starttime
%     Title           - title string as a cell array (def={''})

% set defaults
p.FrameBaseName='frame';
p.ScriptToAdd='none';
p.ImageResolution='-r200';
p.AxisLims=[];  
p.Title='';     
p.ColorMin=NaN;  
p.ColorMax=NaN;  
p.ColorMap=jet(32);  
p.SetUpAxes=true;
p.AxesHandle=NaN;
p.FileName='maxele.63.nc';
p.NewFig=false;
p.Grid='';

echo off
AdcUtil=adcirc_util_init;

p=parse_pv_pairs(p,varargin);

if nargin==0
%   disp('Plot63nc(g,nc) OR:')
%   disp('Plot63nc(fgs,nc,p1,v1,p2,v2,...)')
%   return
    if ~exist(p.FileName,'file')
        error('maxele.63.nc file not found.')
    end   
end

try 
    nc=ncgeodataset(p.FileName);
catch
    error('nc file %s not found or failed to open.',p.FileName);
end

% get filename in case its not maxele
[~,b,~]=fileparts(p.FileName); 
fvar=strtok(b,'.');

if isempty(p.Grid)
    g=ExtractGrid(nc);
else
    g=p.Grid;
end

if p.NewFig
    figure
end

if exist(p.ScriptToAdd,'file')
   eval(p.ScriptToAdd)
end

if ~(isnan(p.ColorMin) && isnan(p.ColorMax))    
    set(gca,'CLim',[p.ColorMin p.ColorMax])
end
colormap(p.ColorMap)
colorbar
   
if isempty(p.Title)
    p.Title=sprintf('%s',p.FileName);
end

zz=nc{AdcUtil.MaxVarMapping.(fvar)}(:)';
zz(zz<AdcUtil.MaxVarNaN.(fvar))=NaN;
h=colormesh2d(g,zz);
hz0=lcontour(g,'z',0,'Color','k');
if ~(isnan(p.ColorMin) && isnan(p.ColorMax))    
   set(gca,'CLim',[p.ColorMin p.ColorMax])
end
plotbnd(g)

%set(gca,'CLim',[-1 1]*max(zz))
%shading flat
%set_height(h,1)
%[hst,axx]=stamp_right(datestr(t(i)));
axis('equal')
%axis('tight')
if ~isnan(p.AxisLims)
   axis(p.AxisLims) 
end

title(p.Title,'FontSize',14,'Interpreter','none')
drawnow
