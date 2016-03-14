function OwiStruct=LoadOwi(FileNamePrefix,Stride,IterEnd,NoHeader)
% OwiStruct=LoadOwi(FileNamePrefix,Stride,IterEnd,NoHeader)

if ~exist('FileNamePrefix','var')
    FileNamePrefix='fort';
    BasinPreFileName=[FileNamePrefix '.221'];
    BasinWinFileName=[FileNamePrefix '.222'];
    RegionPreFileName=[FileNamePrefix '.223'];
    RegionWinFileName=[FileNamePrefix '.224'];
end
if iscell(FileNamePrefix)
    if length(FileNamePrefix) ~=4
        error ('If first arg to LoadOwi is a cell, then it must contain 4 values.')
    end
   BasinPreFileName=FileNamePrefix{1};
   BasinWinFileName=FileNamePrefix{2};
   RegionPreFileName=FileNamePrefix{3};
   RegionWinFileName=FileNamePrefix{4};
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

OwiStruct.Basin=[];
OwiStruct.Region=[];

% Basin Pressure Read
if ~isempty(BasinPreFileName) && exist(BasinPreFileName,'file')
   temp=LoadOwi0(BasinPreFileName,Stride,IterEnd,NoHeader);
   OwiStruct.Basin.time=[temp.time];
   OwiStruct.Basin.iLat=[temp.iLat];
   OwiStruct.Basin.iLong=[temp.iLong];
   OwiStruct.Basin.DX=[temp.DX];
   OwiStruct.Basin.DY=[temp.DY];
   OwiStruct.Basin.SWLat=[temp.SWLat];
   OwiStruct.Basin.SWLon=[temp.SWLon];
   for i=1:length(temp)
      OwiStruct.Basin.XGrid{i}=temp(i).SWLon+(0:temp(i).iLong-1)*temp(i).DX;
      OwiStruct.Basin.YGrid{i}=temp(i).SWLat+(0:temp(i).iLat-1)'*temp(i).DY;
      OwiStruct.Basin.Pre{i}=temp(i).Pre';
   end
else
    if isempty(BasinPreFileName) 
        disp(['No Basin Pre file specified.'])
    else
        disp(['Basin Pre file (' BasinPreFileName ') DNE.'])
    end
end

% Basin Wind Read
if ~isempty(BasinWinFileName) && exist(BasinWinFileName,'file')
   temp=LoadOwi0(BasinWinFileName,Stride,IterEnd,NoHeader);
   OwiStruct.Basin.time=[temp.time];
   OwiStruct.Basin.iLat=[temp.iLat];
   OwiStruct.Basin.iLong=[temp.iLong];
   OwiStruct.Basin.DX=[temp.DX];
   OwiStruct.Basin.DY=[temp.DY];
   OwiStruct.Basin.SWLat=[temp.SWLat];
   OwiStruct.Basin.SWLon=[temp.SWLon];
   maxspd=NaN*ones(size(temp(1).Win.u));
   for i=1:length(temp)
      u=temp(i).Win.u';
      v=temp(i).Win.v';
      OwiStruct.Basin.WinU{i}=u;
      OwiStruct.Basin.WinV{i}=v;
      spd=abs(u+sqrt(-1)*v);
    %  maxspd=max(maxspd,spd);
   end
   if ~isfield(OwiStruct.Basin,'XGrid')
       for i=1:length(temp)
           OwiStruct.Basin.XGrid{i}=temp(i).SWLon+(0:temp(i).iLong-1)*temp(i).DX;
           OwiStruct.Basin.YGrid{i}=temp(i).SWLat+(0:temp(i).iLat-1)'*temp(i).DY;
       end
   end
   OwiStruct.Basin.maxspd=maxspd;
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
   OwiStruct.Region.time=[temp.time];
   OwiStruct.Region.iLat=[temp.iLat];
   OwiStruct.Region.iLong=[temp.iLong];
   OwiStruct.Region.DX=[temp.DX];
   OwiStruct.Region.DY=[temp.DY];
   OwiStruct.Region.SWLat=[temp.SWLat];
   OwiStruct.Region.SWLon=[temp.SWLon];
   for i=1:length(temp)
      OwiStruct.Region.XGrid{i}=temp(i).SWLon+(0:temp(i).iLong-1)*temp(i).DX;
      OwiStruct.Region.YGrid{i}=temp(i).SWLat+(0:temp(i).iLat-1)'*temp(i).DY;
      OwiStruct.Region.Pre{i}=temp(i).Pre';
   end
   FillRegionGrid=false;
else
    if isempty(RegionPreFileName)
        fprintf('No Region Pre file specified.\n')
    else
       disp(['Region Pre file (' RegionPreFileName ') DNE.'])
    end
end

if ~isempty(RegionWinFileName) && exist(RegionWinFileName,'file')
   temp=LoadOwi0(RegionWinFileName,Stride,IterEnd,NoHeader);
   OwiStruct.Region.time=[temp.time];
   OwiStruct.Region.iLat=[temp.iLat];
   OwiStruct.Region.iLong=[temp.iLong];
   OwiStruct.Region.DX=[temp.DX];
   OwiStruct.Region.DY=[temp.DY];
   OwiStruct.Region.SWLat=[temp.SWLat];
   OwiStruct.Region.SWLon=[temp.SWLon];
   maxspd=NaN*ones(size(temp(1).Win.u));
   for i=1:length(temp)
      u=temp(i).Win.u';
      v=temp(i).Win.v';
      OwiStruct.Region.WinU{i}=temp(i).Win.u';
      OwiStruct.Region.WinV{i}=temp(i).Win.v';
      spd=abs(u+sqrt(-1)*v);
      maxspd=max(maxspd,spd');
      
      if FillRegionGrid
        OwiStruct.Region.XGrid{i}=temp(i).SWLon+(0:temp(i).iLong-1)*temp(i).DX;
        OwiStruct.Region.YGrid{i}=temp(i).SWLat+(0:temp(i).iLat-1)'*temp(i).DY;
      end
      
   end
   OwiStruct.Region.maxspd=maxspd;
else
    if isempty(RegionWinFileName)
        fprintf('No Region Win file specified.\n')
    else
       disp(['Region Win file (' RegionWinFileName ') DNE.'])    
    end
end



function D=LoadOwi0(FileName,Stride,IterEnd,NoHeader)
%LOADOWI0
% D=LoadOwi0(FileName,Stride);

Debug=0;

if ~exist(FileName,'file')
   error('%s does not exist.',FileName)
end

[~,~,Type]=fileparts(FileName);
Type=Type(2:end);
if ~any(strcmpi(Type,{'221','222','223','224','pre','wnd','win'}))
   error('WND/PRE file extension (%s) must be 221|222|223|224|win|wnd|pre.',Type)
end
fid=fopen(FileName,'r');

if ~exist('Stride')
   Stride=1;
end

if ~NoHeader
headerline=fgets(fid);
%disp(headerline)
%temp=strread(headerline,'%s');
%[a,b,c,d,e]=strread(temp{4},'%4d%2d%2d%2d%2d');
end

time_string='iLat=%diLong=%dDX=%fDY=%fSWLat=%fSWLon=%fDT=%4d%2d%2d%2d%2d';

k=0;j=0;

if Debug>-1
   fprintf('Reading from %s \n',FileName)
end

while ~feof(fid)
   k=k+1;
   l=fgets(fid);
   %fgetl(fid);
   %[iLat,iLong,DX,DY,SWLat,SWLon,y,m,d,h,mn]
   A=sscanf(l,time_string);
   %disp(A);
   ctime=datenum(A(7),A(8),A(9),A(10),A(11),0);
   iLat=A(1);
   iLong=A(2);
   
   if Debug>0
      fprintf('Scanning [%d x %d] "%s" values at time=%s\n',iLong,iLat,Type,datestr(ctime,0))
   end
   
   u=fscanf(fid,'%f',[iLong,iLat]);  
   
   if any(strcmpi(Type,{'222','224','wnd','win'}))
     if Debug>0
        fprintf('Scanning [%d x %d] "%s" values at time=%s\n',iLong,iLat,Type,datestr(ctime,0))
     end
     v=fscanf(fid,'%f',[iLong,iLat]);
   end
   
   if mod(k-1,Stride)==0
      if Debug>1
         fprintf('   Storing [%d x %d] "%s" values at time=%s\n',iLong,iLat,Type,datestr(ctime,0))
      end
      j=j+1;
      D(j).time=ctime;
      D(j).iLat=iLat;
      D(j).iLong=iLong;
      D(j).DX=A(3);
      D(j).DY=A(4);
      D(j).SWLat=A(5);
      D(j).SWLon=A(6);
      if any(strcmp(Type,{'221','223','pre'}))
         D(j).Pre=u;
      else
         D(j).Win.u=u;
      end
      
      if Debug>1
         fprintf('   Min,Max u = %f %f\n',min(u(:)),max(u(:)))
      end
      if any(strcmpi(Type,{'222','224','wnd','win'}))
         D(j).Win.v=v;
         if Debug>1
            fprintf('   Min,Max v = %f %f\n',min(u(:)),max(u(:)))
         end
      end      
   end
   fgetl(fid);
   if k==IterEnd, break,end
end    
   
fclose(fid);
