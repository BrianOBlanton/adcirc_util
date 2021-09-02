function AnimOwi(OwiStruct,varargin)
% AnimOwi generate an animatiom of an OWI set of files
%  AnimOwi(OwiStruct,varargin)
% Required inputs: 
%     OwiStruct - structure read in by LoadOwi (required)
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
%     BasinVectorColor - def='k';
%     RegionVectorColor - def='w';
%     BasinVectorStride - def=5;
%     RegionVectorStride - def=20;

if nargin==0
   disp('AnimOwi(OwiStruct)  OR: ')
   disp('AnimOwi(OwiStruct,p1,v1,p2,v2,...)')
   return
end

FrameBaseName='frame';
ScriptToAdd='none';
ImageResolution='-r150';
AxisLims=[];  
Title={''};     
IterStart=1; 
IterStride=1;
IterStop=-1;  
ColorMin=935;  
ColorMax=1015;  
ColorMap=jet(16);  

BasinVectorColor='k';
RegionVectorColor='k';

BasinVectorStride=10;
RegionVectorStride=10;

BasinVectorScaleFac=100;
RegionVectorScaleFac=50;

BasinPressureDraw=true;
BasinVectorDraw=true;
RegionPressureDraw=false;
RegionVectorDraw=false;

if ~isfield(OwiStruct.Basin,'Pre')
    BasinPressureDraw=false;
end
if ~isfield(OwiStruct.Basin,'WinU')
    BasinVectorDraw=false;
end
if ~isfield(OwiStruct.Region,'Pre')
    RegionPressureDraw=false;
end
if ~isfield(OwiStruct.Region,'WinU')
    RegionVectorDraw=false;
end


% Strip off propertyname/value pairs in varargin not related to
% "line" object properties.
k=1;
while k<length(varargin),
  switch lower(varargin{k}),
    case 'basinvectorscalefac',
      BasinVectorScaleFac=varargin{k+1};
      varargin([k k+1])=[];
   case 'regionvectorscalefac',
      RegionVectorScaleFac=varargin{k+1};
      varargin([k k+1])=[];
    case 'basinvectorcolor',
      BasinVectorColor=varargin{k+1};
      varargin([k k+1])=[];
   case 'regionvectorcolor',
      RegionVectorColor=varargin{k+1};
      varargin([k k+1])=[];
    case 'basinvectorstride',
      BasinVectorStride=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionvectorstride',
      RegionVectorStride=varargin{k+1};
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
    case 'title',
      Title=varargin{k+1};
      varargin([k k+1])=[];
    case 'iterstart',
      IterStart=varargin{k+1};
      varargin([k k+1])=[];
    case 'iterstride',
      IterStride=varargin{k+1};
      varargin([k k+1])=[];
    case 'iterstop',
      IterStop=varargin{k+1};
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
    otherwise
      k=k+2;
  end;
end;

if length(varargin)<2
   varargin={};
end

Basin=isfield(OwiStruct,'Basin');
%Region=isfield(OwiStruct,'Region');
Region=~isempty(OwiStruct.Region);

time=[OwiStruct.Basin.time];
nTimes=length(time);

if IterStop==-1,IterStop=nTimes;end

if (IterStart<1 || IterStart>nTimes)
   error('IterStart must be between %d and %d',1,nTimes)
elseif (IterStop<1 || IterStop>nTimes)
   error('IterStop must be between %d and %d',1,nTimes)
elseif (IterStop<IterStart)
       error('IterStop must be greater than or equal to IterStart')
end

if IterStride<1
    error('Stride must be greater than 1.')
end

disp(sprintf('Starting iteration =  %d',IterStart))
disp(sprintf('Stride =  %d',IterStride))
disp(sprintf('Stopping iteration =  %d',IterStop))

if ~iscell(Title)
   error('Title to AnimOwi must be a cell array')
end
tlen=length(Title);

% generate initial grid
i=1;
SWLon=OwiStruct.Basin.SWLon(i);
SWLat=OwiStruct.Basin.SWLat(i);
iLong=OwiStruct.Basin.iLong(i);
iLat=OwiStruct.Basin.iLat(i);
DX=OwiStruct.Basin.DX(i);
DY=OwiStruct.Basin.DY(i);
x=SWLon+(0:iLong-1)*DX;
y=SWLat+(0:iLat-1)*DY;
[Xb,Yb]=meshgrid(x,y);

if isempty(AxisLims)
    minX=min(x);maxX=max(x);
    minY=min(y);maxY=max(y);
else
    minX=AxisLims(1);
    maxX=AxisLims(2);
    minY=AxisLims(3);
    maxY=AxisLims(4);
end

% set up figure
if (1)
   figure
   line(Xb(1,:),Yb(1,:));
   line(Xb(end,:),Yb(end,:));
   line(Xb(:,1),Yb(:,1));
   line(Xb(:,end),Yb(:,end));
   plotcoast('states')
   plotcoast('worldcoast')
   set(gca,'FontSize',14)
end
axis('equal')
axis([minX maxX minY maxY])
%plot_google_map('Maptype','hybrid') % ,'scale',2)

if exist(ScriptToAdd,'file')
   eval(ScriptToAdd)
end

if BasinPressureDraw || RegionPressureDraw
    set(gca,'CLim',[ColorMin ColorMax])
    colormap(ColorMap)
    hcb=colorbar;
    mainax=gca;
    hcb.Label.String='Atmospheric Pressure [millibars]';
    %axes(mainax)
end

set(gcf,'Renderer','opengl')

% vidObj=VideoWriter('hbl.anim.avi');
% vidObj.FrameRate = 6;
% vidObj.Quality=50;
% open(vidObj);

hold on

for i=IterStart:IterStride:IterStop
   
   hpb=[];
   hvb=[];
   hpr=[];
   hvr=[];
   hlr=[];
   hc=[];
   if Basin
      
      SWLon=OwiStruct.Basin.SWLon(i);
      SWLat=OwiStruct.Basin.SWLat(i);
      iLong=OwiStruct.Basin.iLong(i);
      iLat=OwiStruct.Basin.iLat(i);
      DX=OwiStruct.Basin.DX(i);
      DY=OwiStruct.Basin.DY(i);
      x=SWLon+(0:iLong-1)*DX;
      y=SWLat+(0:iLat-1)*DY;
      [Xb,Yb]=meshgrid(x,y);
      
      if BasinPressureDraw
          p=OwiStruct.Basin.Pre{i};
          hpb=pcolor(Xb,Yb,p);
          %set(hpb,'FaceAlpha',.5) 
          shading flat
      end
      
      u=OwiStruct.Basin.WinU{i};
      v=OwiStruct.Basin.WinV{i};
      [s,d]=vecmagdir(u,v);
      hvb=vecplot(Xb,Yb,u,v,'ScaleLabel','no scale',...
          'Stride',BasinVectorStride,...
          'ScaleFac',BasinVectorScaleFac,...
          'Color',BasinVectorColor);
      set_height(hvb,10)
      %[c,hc]=contour(Xb,Yb,s*2.23,[25:25:150]);
      %set(hc,'Color','g','Linewidth',2) 
   
   end   % end if Basin
   
%    if Region
%       
%       SWLon=OwiStruct.Region.SWLon(i);
%       SWLat=OwiStruct.Region.SWLat(i);
%       iLong=OwiStruct.Region.iLong(i);
%       iLat=OwiStruct.Region.iLat(i);
%       DX=OwiStruct.Region.DX(i);
%       DY=OwiStruct.Region.DY(i);
%       x=SWLon+(0:iLong-1)*DX;
%       y=SWLat+(0:iLat-1)*DY;
%       [Xr,Yr]=meshgrid(x,y);
%       u=OwiStruct.Region.WinU{i};
%       v=OwiStruct.Region.WinV{i};
%       
%       hlr(1)=line(Xr(1,:),Yr(1,:));
%       hlr(2)=line(Xr(end,:),Yr(end,:));
%       hlr(3)=line(Xr(:,1),Yr(:,1));
%       hlr(4)=line(Xr(:,end),Yr(:,end));
%    
%       if RegionPressureDraw
%         p=OwiStruct.Region.Pre{i};
%         hpr=pcolor(Xr,Yr,p);
%         set(gco,'FaceAlpha',.5) 
%         %shading flat
%       end
%       
%       hvr=vecplot(Xr,Yr,u,v,'ScaleLabel','no scale',...
%           'Stride',RegionVectorStride,...
%           'ScaleFac',RegionVectorScaleFac,...
%           'Color',RegionVectorColor);
%       
%    end   % end if Region
%    
%   Title{tlen+1}=sprintf('Day %s',sprintf('%4.2f',(time(i)-datenum(2008,8,28))));
   Title{tlen+1}=datestr(time(i),0);
   title(Title,'FontSize',18);
   
   drawnow
   
   %pause
   
%   currFrame=getframe(gcf);
%   writeVideo(vidObj,currFrame);
   
   fnamebase=sprintf('%s_%03d',FrameBaseName,i);
   disp(sprintf('Printing %s',fnamebase))
   print('-dpng',ImageResolution,sprintf('%s.png',fnamebase));
   ConvertImage(fnamebase,'png','gif')
   
   delete(hpb)
   delete(hvb)
   delete(hpr)
   delete(hvr)
   delete(hlr)
   delete(hc)

end

close(vidObj);

