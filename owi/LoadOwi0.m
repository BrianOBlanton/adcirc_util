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
