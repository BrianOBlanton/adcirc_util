function [h,g]=PlotMax63nc(varargin)
% [h,g]=PlotMax63nc(varargin)
%
% Plot63nc draw a timestep of an ADCIRC 63 file in netCDF format
% Required inputs: 
% 
% P/V pairs:
%     AxisLims        - plot axis limits, as in the axis command (def=grid lims)
%     Title           - title string as a cell array (def={''})
%     ColorMin        - min pressure to clip below (def=950)
%     ColorMax        - max pressure to clip above (def=1030)
%     ColorMap        - colormap to use (def=jet(32))
%     ScriptToAdd     - script that defines plot overlays (def='none')
%     ImageResolution - (def='-r200';)
%     FrameBaseName   - base of image output file name (def='frame')
%     StartingTime    - datenum of starttime

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
p.filename='maxele.63.nc';
p.NewFig=false;
p.Grid='';

echo off
AdcUtil=adcirc_util_init;

p=parse_pv_pairs(p,varargin);

fvar=strtok(p.filename,'.');

if nargin==0
%   disp('Plot63nc(g,nc) OR:')
%   disp('Plot63nc(fgs,nc,p1,v1,p2,v2,...)')
%   return
    if ~exist(p.filename,'file')
        error('If no arguments, maxele.63.nc must locally exist.')
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
    p.Title=sprintf('%s',p.filename);
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
%axis('equal')
%axis('tight')
title(p.Title,'FontSize',14)
drawnow
