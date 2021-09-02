function owi=LoadOwi(FileNamePrefix,Stride,IterEnd,NoHeader)
% Load OWI formatted wind/pressure files into MATLAB
%
% owi=LoadOwi(FileNamePrefix,Stride,IterEnd,NoHeader)
%
% FileNamePrefix - filename if other than "fort".
% if each file is named differently, pass a cell array of filenames.
%

if ~exist('FileNamePrefix','var'),FileNamePrefix='fort';end

if iscell(FileNamePrefix)
    if length(FileNamePrefix) ~=6
        error ('If first arg to LoadOwi is a cell, then it must contain 6 values.')
    end
    BasinPreFileName=FileNamePrefix{1};
    BasinWinFileName=FileNamePrefix{2};
    RegionPreFileName=FileNamePrefix{3};
    RegionWinFileName=FileNamePrefix{4};
    LocalPreFileName=FileNamePrefix{3};
    LocalWinFileName=FileNamePrefix{4};
else
    BasinPreFileName=[FileNamePrefix '.221'];
    BasinWinFileName=[FileNamePrefix '.222'];
    RegionPreFileName=[FileNamePrefix '.223'];
    RegionWinFileName=[FileNamePrefix '.224'];
    LocalPreFileName=[FileNamePrefix '.225'];
    LocalWinFileName=[FileNamePrefix '.226'];
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
owi.Local=[];

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
      owi.Basin.WinSpd{i}=spd;
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


% Region Pressure Read
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

% Region Wind Read
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
      owi.Region.WinSpd{i}=spd;

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


% Local Pressure Read
FillLocalGrid=true;
if ~isempty(LocalPreFileName) && exist(LocalPreFileName,'file')
   temp=LoadOwi0(LocalPreFileName,Stride,IterEnd,NoHeader);
   owi.Local.time=[temp.time];
   owi.Local.iLat=[temp.iLat];
   owi.Local.iLong=[temp.iLong];
   owi.Local.DX=[temp.DX];
   owi.Local.DY=[temp.DY];
   owi.Local.SWLat=[temp.SWLat];
   owi.Local.SWLon=[temp.SWLon];
   minp=NaN*ones(size(temp(1).Pre));
   for i=1:length(temp)
      owi.Local.XGrid{i}=temp(i).SWLon+(0:temp(i).iLong-1)*temp(i).DX;
      owi.Local.YGrid{i}=temp(i).SWLat+(0:temp(i).iLat-1)'*temp(i).DY;
      owi.Local.Pre{i}=temp(i).Pre';
      minp=min(minp,owi.Local.Pre{i}');
   end
   FillLocalGrid=false;
   owi.Local.minp=minp;
else
    if isempty(LocalPreFileName)
        fprintf('No Local Pre file specified.\n')
    else
       disp(['Local Pre file (' LocalPreFileName ') DNE.'])
    end
end

% Local Wind Read
if ~isempty(LocalWinFileName) && exist(LocalWinFileName,'file')
   temp=LoadOwi0(LocalWinFileName,Stride,IterEnd,NoHeader);
   owi.Local.time=[temp.time];
   owi.Local.iLat=[temp.iLat];
   owi.Local.iLong=[temp.iLong];
   owi.Local.DX=[temp.DX];
   owi.Local.DY=[temp.DY];
   owi.Local.SWLat=[temp.SWLat];
   owi.Local.SWLon=[temp.SWLon];
   maxspd=NaN*ones(size(temp(1).Win.u));
   for i=1:length(temp)
      u=temp(i).Win.u';
      v=temp(i).Win.v';
      owi.Local.WinU{i}=temp(i).Win.u';
      owi.Local.WinV{i}=temp(i).Win.v';
      spd=abs(u+sqrt(-1)*v);
      owi.Local.WinSpd{i}=spd;

      maxspd=max(maxspd,spd');
      
      if FillLocalGrid
        owi.Local.XGrid{i}=temp(i).SWLon+(0:temp(i).iLong-1)*temp(i).DX;
        owi.Local.YGrid{i}=temp(i).SWLat+(0:temp(i).iLat-1)'*temp(i).DY;
      end
      
   end
   owi.Local.maxspd=maxspd;
else
    if isempty(LocalWinFileName)
        fprintf('No Local Win file specified.\n')
    else
       disp(['Local Win file (' LocalWinFileName ') DNE.'])    
    end
end

