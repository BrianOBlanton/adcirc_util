function owi=LoadOwi(FileNamePrefix,Stride,IterEnd,NoHeader)
% Load OWI formatted wind/pressure files into MATLAB
%
% owi=LoadOwi(FileNamePrefix,Stride,IterEnd,NoHeader)
%
% FileNamePrefix - filename if other than "fort".  

if ~exist('FileNamePrefix'),FileNamePrefix='fort';end

if iscell(FileNamePrefix)
    if length(FileNamePrefix) ~=4
        error ('If first arg to LoadOwi is a cell, then it must contain 4 values.')
    end
    BasinPreFileName=FileNamePrefix{1};
    BasinWinFileName=FileNamePrefix{2};
    RegionPreFileName=FileNamePrefix{3};
    RegionWinFileName=FileNamePrefix{4};
else
    BasinPreFileName=[FileNamePrefix '.221'];
    BasinWinFileName=[FileNamePrefix '.222'];
    RegionPreFileName=[FileNamePrefix '.223'];
    RegionWinFileName=[FileNamePrefix '.224'];
end

if ~exist('Stride','var')
    Stride=1;
end
if ~exist('IterEnd','var')
    IterEnd=-1;
end
if ~exist('NoHeader','var')
    NoHeader=false;
end

owi.Basin=[];
owi.Region=[];

% Basin Pressure Read
if ~isempty(BasinPreFileName) && exist(BasinPreFileName,'file')
   temp=LoadOwi0(BasinPreFileName,Stride,IterEnd,NoHeader);
   owi.Basin.time=[temp.time];
   owi.Basin.iLat=[temp.iLat];
   owi.Basin.iLong=[temp.iLong];
   owi.Basin.DX=[temp.DX];
   owi.Basin.DY=[temp.DY];
   owi.Basin.SWLat=[temp.SWLat];
   owi.Basin.SWLon=[temp.SWLon];
   minp=NaN*ones(size(temp(1).Pre));
   for i=1:length(temp)
      owi.Basin.XGrid{i}=temp(i).SWLon+(0:temp(i).iLong-1)*temp(i).DX;
      owi.Basin.YGrid{i}=temp(i).SWLat+(0:temp(i).iLat-1)'*temp(i).DY;
      owi.Basin.Pre{i}=temp(i).Pre';
      minp=min(minp,owi.Basin.Pre{i}');
   end
   owi.Basin.minp=minp;
else
    if isempty(BasinPreFileName) 
        disp('No Basin Pre file specified.')
    else
        disp(['Basin Pre file (' BasinPreFileName ') DNE.'])
    end
end

% Basin Wind Read
if ~isempty(BasinWinFileName) && exist(BasinWinFileName,'file')
   temp=LoadOwi0(BasinWinFileName,Stride,IterEnd,NoHeader);
   owi.Basin.time=[temp.time];
   owi.Basin.iLat=[temp.iLat];
   owi.Basin.iLong=[temp.iLong];
   owi.Basin.DX=[temp.DX];
   owi.Basin.DY=[temp.DY];
   owi.Basin.SWLat=[temp.SWLat];
   owi.Basin.SWLon=[temp.SWLon];
   maxspd=NaN*ones(size(temp(1).Win.u));
   for i=1:length(temp)
      u=temp(i).Win.u';
      v=temp(i).Win.v';
      owi.Basin.WinU{i}=u;
      owi.Basin.WinV{i}=v;
      spd=abs(u+sqrt(-1)*v);
      maxspd=max(maxspd,spd');
   end
   if ~isfield(owi.Basin,'XGrid')
       for i=1:length(temp)
           owi.Basin.XGrid{i}=temp(i).SWLon+(0:temp(i).iLong-1)*temp(i).DX;
           owi.Basin.YGrid{i}=temp(i).SWLat+(0:temp(i).iLat-1)'*temp(i).DY;
       end
   end
   owi.Basin.maxspd=maxspd;
else
    if isempty(BasinWinFileName)
        fprintf('No Basin Win file specified.\n')
    else
        disp(['Basin Win file (' BasinWinFileName ') DNE.'])
    end
end


FillRegionGrid=true;
if ~isempty(RegionPreFileName) && exist(RegionPreFileName,'file')
   temp=LoadOwi0(RegionPreFileName,Stride,IterEnd,NoHeader);
   owi.Region.time=[temp.time];
   owi.Region.iLat=[temp.iLat];
   owi.Region.iLong=[temp.iLong];
   owi.Region.DX=[temp.DX];
   owi.Region.DY=[temp.DY];
   owi.Region.SWLat=[temp.SWLat];
   owi.Region.SWLon=[temp.SWLon];
   minp=NaN*ones(size(temp(1).Pre));
   for i=1:length(temp)
      owi.Region.XGrid{i}=temp(i).SWLon+(0:temp(i).iLong-1)*temp(i).DX;
      owi.Region.YGrid{i}=temp(i).SWLat+(0:temp(i).iLat-1)'*temp(i).DY;
      owi.Region.Pre{i}=temp(i).Pre';
      minp=min(minp,owi.Region.Pre{i}');
   end
   FillRegionGrid=false;
   owi.Region.minp=minp;
else
    if isempty(RegionPreFileName)
        fprintf('No Region Pre file specified.\n')
    else
       disp(['Region Pre file (' RegionPreFileName ') DNE.'])
    end
end

if ~isempty(RegionWinFileName) && exist(RegionWinFileName,'file')
   temp=LoadOwi0(RegionWinFileName,Stride,IterEnd,NoHeader);
   owi.Region.time=[temp.time];
   owi.Region.iLat=[temp.iLat];
   owi.Region.iLong=[temp.iLong];
   owi.Region.DX=[temp.DX];
   owi.Region.DY=[temp.DY];
   owi.Region.SWLat=[temp.SWLat];
   owi.Region.SWLon=[temp.SWLon];
   maxspd=NaN*ones(size(temp(1).Win.u));
   for i=1:length(temp)
      u=temp(i).Win.u';
      v=temp(i).Win.v';
      owi.Region.WinU{i}=temp(i).Win.u';
      owi.Region.WinV{i}=temp(i).Win.v';
      spd=abs(u+sqrt(-1)*v);
      maxspd=max(maxspd,spd');
      
      if FillRegionGrid
        owi.Region.XGrid{i}=temp(i).SWLon+(0:temp(i).iLong-1)*temp(i).DX;
        owi.Region.YGrid{i}=temp(i).SWLat+(0:temp(i).iLat-1)'*temp(i).DY;
      end
      
   end
   owi.Region.maxspd=maxspd;
else
    if isempty(RegionWinFileName)
        fprintf('No Region Win file specified.\n')
    else
       disp(['Region Win file (' RegionWinFileName ') DNE.'])    
    end
end

