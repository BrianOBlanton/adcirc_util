function h=PlotOwi(OwiStruct,i,varargin)
%
% h=PlotOwi(OwiStruct,index,p1,v1,...)
%
% OwiStruct - reguired, as loaded using LoadOwi
% index - required, time index to draw
%
% Optional PV args/defaults:
%
%     DrawScale=false;
%
%     ContourSpeed=false;
%     
%     ColorMin=980;  
%     ColorMax=1030;  
%     ColorMap=jet(32);  
%     
%     BasinVectorColor='r';
%     RegionVectorColor='b';
%     LocalVectorColor='g';
%   
%     BasinVectorStride=10;
%     RegionVectorStride=10;
%     LocalVectorStride=10;
%     
%     BasinVectorScaleFac=50;
%     RegionVectorScaleFac=50;
%     LocalVectorScaleFac=50;
%     
%     BasinPressureDraw=true;
%     BasinVectorDraw=true;
%     RegionPressureDraw=true;
%     RegionVectorDraw=true;
%     LocalPressureDraw=true;
%     LocalVectorDraw=true;
%

% Default propertyname values

DrawScale=false;

ContourSpeed=false;
ColorMin=980;  
ColorMax=1030;  
ColorMap=jet(32);  

BasinVectorColor='k';
RegionVectorColor='k';
LocalVectorColor='k';

BasinVectorStride=10;
RegionVectorStride=10;
LocalVectorStride=10;

BasinVectorScaleFac=50;
RegionVectorScaleFac=50;
LocalVectorScaleFac=50;

BasinVectorScaleLabel='m/s';
RegionVectorScaleLabel='m/s';
LocalVectorScaleLabel='m/s';

BasinPressureDraw=true;
BasinVectorDraw=true;
RegionPressureDraw=true;
RegionVectorDraw=true;
LocalPressureDraw=true;
LocalVectorDraw=true;

ScaleXor=[];
ScaleYor=[];

AxisLims=[];

% Strip off propertyname/value pairs in varargin not related to
% "line" object properties.
k=1;
while k<length(varargin)
  switch lower(varargin{k})
    case 'contourspeed'
      ContourSpeed=varargin{k+1};
      varargin([k k+1])=[];    
      ColorMin=0;
      ColorMax=50;
    case 'basinpressuredraw'
      BasinPressureDraw=varargin{k+1};
      varargin([k k+1])=[];
    case 'basinvectorscalefac'
      BasinVectorScaleFac=varargin{k+1};
      varargin([k k+1])=[];
    case 'basinvectordraw'
      BasinVectorDraw=varargin{k+1};
      varargin([k k+1])=[];
    case 'basinvectorstride'
      BasinVectorStride=varargin{k+1};
      varargin([k k+1])=[];
    case 'basinvectorcolor'
      BasinVectorColor=varargin{k+1};
      varargin([k k+1])=[];
    case 'basinvectorscalelabel'
      BasinVectorScaleLabel=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionpressuredraw'
      RegionPressureDraw=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionvectordraw'
      RegionVectorDraw=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionvectorstride'
      RegionVectorStride=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionvectorscalefac'
      RegionVectorScaleFac=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionvectorcolor'
      RegionVectorColor=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionvectorscalelabel'
      RegionVectorScaleLabel=varargin{k+1};
      varargin([k k+1])=[];
    case 'localpressuredraw'
      LocalPressureDraw=varargin{k+1};
      varargin([k k+1])=[];
    case 'localvectordraw'
      LocalVectorDraw=varargin{k+1};
      varargin([k k+1])=[];
    case 'localvectorstride'
      LocalVectorStride=varargin{k+1};
      varargin([k k+1])=[];
    case 'localvectorscalefac'
      LocalVectorScaleFac=varargin{k+1};
      varargin([k k+1])=[];
    case 'localvectorscalelabel'
      LocalVectorScaleLabel=varargin{k+1};
      varargin([k k+1])=[];
    case 'localvectorcolor'
      LocalVectorColor=varargin{k+1};
      varargin([k k+1])=[];
    case 'axislims'
      AxisLims=varargin{k+1};
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
    case 'scalexor'
      ScaleXor=varargin{k+1};
      varargin([k k+1])=[];
    case 'scaleyor'
      ScaleYor=varargin{k+1};
      varargin([k k+1])=[];
    case 'drawscale'
      DrawScale=varargin{k+1};
      varargin([k k+1])=[];
    otherwise
      k=k+2;
  end
end

if length(varargin)<2
   varargin={};
end

h=[];

if ContourSpeed
    if isfield(OwiStruct,'Basin')
        if ~isfield(OwiStruct.Basin,'WinSpd')
            error('ContourSpeed is true, but Basin does not contain WinSpd field.')
        end
    elseif isfield(OwiStruct,'Region')
        if ~isfield(OwiStruct.Region,'WinSpd')
            error('ContourSpeed is true, but Region does not contain WinSpd field.')
        end
    elseif isfield(OwiStruct,'Local')
        if ~isfield(OwiStruct.Local,'WinSpd')
            error('ContourSpeed is true, but Local does not contain WinSpd field.')
        end
    end
end

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
if ~isfield(OwiStruct.Local,'Pre')
    LocalPressureDraw=false;
end
if ~isfield(OwiStruct.Local,'WinU')
    LocalVectorDraw=false;
end

% [BasinPressureDraw BasinVectorDraw  RegionPressureDraw RegionVectorDraw]

if isempty(AxisLims)
    if (BasinPressureDraw || BasinVectorDraw)
        SWLon=OwiStruct.Basin.SWLon(1);
        SWLat=OwiStruct.Basin.SWLat(1);
        iLong=OwiStruct.Basin.iLong(1);
        iLat=OwiStruct.Basin.iLat(1);
        DX=OwiStruct.Basin.DX(1);
        DY=OwiStruct.Basin.DY(1);
    else
        SWLon=OwiStruct.Region.SWLon(1);
        SWLat=OwiStruct.Region.SWLat(1);
        iLong=OwiStruct.Region.iLong(1);
        iLat=OwiStruct.Region.iLat(1);
        DX=OwiStruct.Region.DX(1);
        DY=OwiStruct.Region.DY(1);
    end
    x=SWLon+(0:iLong-1)*DX;
    y=SWLat+(0:iLat-1)*DY;
    minX=min(x);maxX=max(x);
    minY=min(y);maxY=max(y);
else
    minX=AxisLims(1);
    maxX=AxisLims(2);
    minY=AxisLims(3);
    maxY=AxisLims(4);
end

%[BasinPressureDraw BasinVectorDraw RegionPressureDraw RegionVectorDraw]
%[minX maxX minY maxY]
bndlines={'LineWidth',1,'Color','k'};
if ~DrawScale && isempty(ScaleXor)
    ScaleAlreadyDrawn=true;
else
    ScaleAlreadyDrawn=false;
end

if (BasinPressureDraw || BasinVectorDraw)
    if isfield(OwiStruct,'Basin')
        
        SWLon=OwiStruct.Basin.SWLon(i);
        SWLat=OwiStruct.Basin.SWLat(i);
        iLong=OwiStruct.Basin.iLong(i);
        iLat=OwiStruct.Basin.iLat(i);
        DX=OwiStruct.Basin.DX(i);
        DY=OwiStruct.Basin.DY(i);
        x=SWLon+(0:iLong-1)*DX;
        y=SWLat+(0:iLat-1)*DY;
        [Xb,Yb]=meshgrid(x,y);
        
        hpb=[];
        hvb=[];
        
        if BasinPressureDraw
            if ContourSpeed
                p=OwiStruct.Basin.WinSpd{i};
            else
                p=OwiStruct.Basin.Pre{i};
            end
            hpb=pcolor(Xb,Yb,p);
            shading flat
            colormap(ColorMap)
            set(gca,'CLim',[ColorMin ColorMax])
        end
        
        axis('equal')
        axis([minX maxX minY maxY])
        
        if BasinVectorDraw
            
            if ScaleAlreadyDrawn
                ScaleXor=[];
                ScaleYor=[];
                BasinVectorScaleLabel='no scale';
            end
            
            
            u=OwiStruct.Basin.WinU{i};
            v=OwiStruct.Basin.WinV{i};
            hvb=vecplot(Xb,Yb,u,v,...
                'Stride',BasinVectorStride,...
                'ScaleFac',BasinVectorScaleFac,...
                'ScaleType','fixed',...
                'ScaleXor',ScaleXor,'ScaleYor',ScaleYor,...
                'Color',BasinVectorColor,...
                'ScaleLabel',BasinVectorScaleLabel);
            ScaleAlreadyDrawn=true;
            
        end
        
        hlb(1)=line(Xb(1,:),Yb(1,:),bndlines{:});
        hlb(2)=line(Xb(end,:),Yb(end,:),bndlines{:});
        hlb(3)=line(Xb(:,1),Yb(:,1),bndlines{:});
        hlb(4)=line(Xb(:,end),Yb(:,end),bndlines{:});
        
        h=[h; hpb; hvb];
        title(string(datetime(datevec(OwiStruct.Basin.time(i)))))
    end   % end if Basin
end

if (RegionPressureDraw || RegionVectorDraw)
    if isfield(OwiStruct,'Region')
        
        SWLon=OwiStruct.Region.SWLon(i);
        SWLat=OwiStruct.Region.SWLat(i);
        iLong=OwiStruct.Region.iLong(i);
        iLat=OwiStruct.Region.iLat(i);
        DX=OwiStruct.Region.DX(i);
        DY=OwiStruct.Region.DY(i);
        x=SWLon+(0:iLong-1)*DX;
        y=SWLat+(0:iLat-1)*DY;
        [Xr,Yr]=meshgrid(x,y);
        u=OwiStruct.Region.WinU{i};
        v=OwiStruct.Region.WinV{i};
        
        hpr=[];hvr=[];
        
        if RegionPressureDraw
            if ~ishold
                hold on
            end
            p=OwiStruct.Region.Pre{i};
            hpr=pcolor(Xr,Yr,p);
            shading flat
            colormap(ColorMap)

            set(gca,'CLim',[ColorMin ColorMax])
        end
        
        if RegionVectorDraw
            if ScaleAlreadyDrawn
                ScaleXor=[];
                ScaleYor=[];
                RegionVectorScaleLabel='no scale';
            end
            hvr=vecplot(Xr,Yr,u,v,...
                'ScaleLabel',RegionVectorScaleLabel,...
                'Stride',RegionVectorStride,...
                'ScaleXor',ScaleXor,'ScaleYor',ScaleYor,...
                'ScaleFac',RegionVectorScaleFac,...
                'Color',RegionVectorColor);
        end
        
        
        hlr(1)=line(Xr(1,:),Yr(1,:),bndlines{:});
        hlr(2)=line(Xr(end,:),Yr(end,:),bndlines{:});
        hlr(3)=line(Xr(:,1),Yr(:,1),bndlines{:});
        hlr(4)=line(Xr(:,end),Yr(:,end),bndlines{:});
        
        h=[h; hpr; hvr];
        
    end   % end if Region
end


if (LocalPressureDraw || LocalVectorDraw)
    if isfield(OwiStruct,'Local')
        
        SWLon=OwiStruct.Local.SWLon(i);
        SWLat=OwiStruct.Local.SWLat(i);
        iLong=OwiStruct.Local.iLong(i);
        iLat=OwiStruct.Local.iLat(i);
        DX=OwiStruct.Local.DX(i);
        DY=OwiStruct.Local.DY(i);
        x=SWLon+(0:iLong-1)*DX;
        y=SWLat+(0:iLat-1)*DY;
        [Xl,Yl]=meshgrid(x,y);
        u=OwiStruct.Local.WinU{i};
        v=OwiStruct.Local.WinV{i};
        
        hpr=[];hvr=[];
        
        if LocalPressureDraw
            if ~ishold
                hold on
            end
            p=OwiStruct.Local.Pre{i};
            hpr=pcolor(Xl,Yl,p);
            shading flat
            colormap(ColorMap)

            set(gca,'CLim',[ColorMin ColorMax])
        end
        
        if LocalVectorDraw
            if ScaleAlreadyDrawn 
                ScaleXor=[];
                ScaleYor=[];
                LocalVectorScaleLabel='no scale';
            end
            hvr=vecplot(Xl,Yl,u,v,...
                'ScaleLabel',LocalVectorScaleLabel,...
                'Stride',LocalVectorStride,...
                'ScaleXor',ScaleXor,'ScaleYor',ScaleYor,...
                'ScaleFac',LocalVectorScaleFac,...
                'Color',LocalVectorColor);
        end
        
        hlr(1)=line(Xl(1,:),Yl(1,:),bndlines{:});
        hlr(2)=line(Xl(end,:),Yl(end,:),bndlines{:});
        hlr(3)=line(Xl(:,1),Yl(:,1),bndlines{:});
        hlr(4)=line(Xl(:,end),Yl(:,end),bndlines{:});
  

        h=[h; hpr; hvr];
        
    end   % end if Local
end
