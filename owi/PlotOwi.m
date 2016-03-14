function h=PlotOwi(OwiStruct,i,varargin)
% h=PlotOwi(OwiStruct,index)

% Default propertyname values


ColorMin=910;  
ColorMax=1030;  
ColorMap=jet(32);  

BasinVectorColor='r';
RegionVectorColor='b';

BasinVectorStride=10;
RegionVectorStride=10;

BasinVectorScaleFac=50;
RegionVectorScaleFac=20;

BasinPressureDraw=true;
BasinVectorDraw=true;
RegionPressureDraw=true;
RegionVectorDraw=true;

AxisLims=[];

% Strip off propertyname/value pairs in varargin not related to
% "line" object properties.
k=1;
while k<length(varargin),
  switch lower(varargin{k}),
    case 'basinpressuredraw',
      BasinPressureDraw=varargin{k+1};
      varargin([k k+1])=[];
    case 'basinvectorscalefac',
      BasinVectorScaleFac=varargin{k+1};
      varargin([k k+1])=[];
    case 'basinvectordraw',
      BasinVectorDraw=varargin{k+1};
      varargin([k k+1])=[];
    case 'basinvectorstride',
      BasinVectorStride=varargin{k+1};
      varargin([k k+1])=[];
    case 'basinvectorcolor',
      BasinVectorColor=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionpressuredraw',
      RegionPressureDraw=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionvectordraw',
      RegionVectorDraw=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionvectorstride',
      RegionVectorStride=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionvectorscalefac',
      RegionVectorScaleFac=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionvectorstride',
      RegionVectorStride=varargin{k+1};
      varargin([k k+1])=[];
    case 'regionvectorcolor',
      RegionVectorColor=varargin{k+1};
      varargin([k k+1])=[];
    case 'axislims',
      AxisLims=varargin{k+1};
      varargin([k k+1])=[];
    case 'colormin',
      ColorMin=varargin{k+1};
      varargin([k k+1])=[];
    case 'colormax',
      ColorMax=varargin{k+1};
      varargin([k k+1])=[];
    otherwise
      k=k+2;
  end;
end;

if length(varargin)<2
   varargin={};
end


k=1;
h=[];

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
        
        hlb(1)=line(Xb(1,:),Yb(1,:));
        hlb(2)=line(Xb(end,:),Yb(end,:));
        hlb(3)=line(Xb(:,1),Yb(:,1));
        hlb(4)=line(Xb(:,end),Yb(:,end));
        
        hpb=[];
        hvb=[];
        
        if BasinPressureDraw
            p=OwiStruct.Basin.Pre{i};
            hpb=pcolor(Xb,Yb,p);
            shading flat
            set(gca,'CLim',[ColorMin ColorMax])
        end
        
       axis([minX maxX minY maxY])
        
        if BasinVectorDraw
            u=OwiStruct.Basin.WinU{i};
            v=OwiStruct.Basin.WinV{i};
            hvb=vecplot(Xb,Yb,u,v,...
                'Stride',BasinVectorStride,...
                'ScaleFac',BasinVectorScaleFac,...
                'Color',BasinVectorColor,...
                'ScaleLabel','no scale');
        end
        
        h=[h; hpb; hvb];
        
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
        
        hlr(1)=line(Xr(1,:),Yr(1,:));
        hlr(2)=line(Xr(end,:),Yr(end,:));
        hlr(3)=line(Xr(:,1),Yr(:,1));
        hlr(4)=line(Xr(:,end),Yr(:,end));
        

        hpr=[];hvr=[];
        
        if RegionPressureDraw
            if ~ishold
                hold on
            end
            p=OwiStruct.Region.Pre{i};
            hpr=pcolor(Xr,Yr,p);
            shading flat
            set(gca,'CLim',[ColorMin ColorMax])
        end
        
        if RegionVectorDraw
            hvr=vecplot(Xr,Yr,u,v,...
                'ScaleLabel','no scale',...
                'Stride',RegionVectorStride,...
                'ScaleFac',RegionVectorScaleFac,...
                'Color',RegionVectorColor);
        end
        h=[h; hpr; hvr];
        
    end   % end if Region
end


