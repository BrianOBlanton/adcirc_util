function WriteOwi(D,filename)
%WriteOwi - output OWI-formatted wind/pressure files from an OwiStruct 
%
% WriteOwi(Owi,filename)
% 
% INPUT:    Owi      - an OwiStruct .  See OwiStruct for details.
%           filename - prefix for file set.  
%                      If empty, WriteOwi writes to fort.22{1,2,3,4} as needed. 
%                      If a string, WriteOwi appends "win", "pre", "_Basin" 
%                      and "_Region" as needed.
%                      To specify the filename for each OWI file, pass in a 
%                      cell with 4 values.  E.g.:
%                      {'hur_Basin.pre','hur_Basin.win','hur_Region.pre','hur_Region.win'}
%

flds={'time','iLat','iLong', 'DX','DY','SWLat','SWLon',...
     'XGrid','YGrid','Pre','WinU', 'WinV'};
 
% check for required fields
if ~isfield(D,'Basin')
    error('OwiStruct must have Basin grid defined.')
else
    if ~all(isfield(D.Basin,flds))
        error('Input OwiStruct.Basin is missing required fields.')
    end
end

Region=true;
if ~isfield(D,'Region') || isempty(D.Region)
    fprintf('No region specified.  So no 223,224 files.\n')
    Region=false;
else
    if ~all(isfield(D.Region,flds))
        error('Input OwiStruct.Region is missing required fields.')
    end
end

if ~exist('filename')
   filename='fort';
   fprintf('No output file name specified, so using "fort".\n')
   % exit if fort.22{1,2,3,4} files already exist.  Will NOT overwrite...
   if  any([exist([filename '.221']) exist([filename '.222'])  ...
            exist([filename '.223']) exist([filename '.224'])])
       str={'Unfortunately at least one fort.22{1,2,3,4} file exists and I cowardly refuse to overwrite.'};
       error(str{:});
   end
   BasinPreFile='fort.221'; 
   BasinWinFile='fort.222';
   RegionPreFile='fort.223'; 
   RegionWinFile='fort.224';
else
    if ischar(filename)
        BasinPreFile=[filename '_Basin.pre'];
        BasinWinFile=[filename '_Basin.win'];
        RegionPreFile=[filename '_Region.pre'];
        RegionWinFile=[filename '_Region.win'];
    elseif iscell(filename)
        if length(filename)==4
            BasinPreFile=filename{1};
            BasinWinFile=filename{2};
            RegionPreFile=filename{3};
            RegionWinFile=filename{4};
        else
            error('Cant make sense of filename input.  Must be length==4 if a cell.')
        end
    end
end
%             1234567890
time_string ='iLat=%4diLong=%4dDX=%6.4fDY=%6.4fSWLat=%8.5fSWLon=%8.4fDT=%12s';
value_string=' %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n';
%value_string=' %9d %9d %9d %9d %9d %9d %9d %9d\n';

bt0=D.Basin.time(1);
bt1=D.Basin.time(end);

% compare times if Region
if Region

    rt0=D.Region.time(1);
    rt1=D.Region.time(end);
    if (bt0 ~= rt0) && (bt1~=rt1)
        error('Time span between Basin and Region are not equal.')
    end
    if length(D.Basin.time) ~= length(D.Region.time)
        error('Number of Time levels between Basin and Region is not equal.')
    end
end

t1=datestr(bt0,30);
t1([9 14 15])=[];
t2=datestr(bt1,30);
t2([9 14 15])=[];
header=sprintf('Oceanweather WIN/PRE Format                        %12s     %12s',t1,t2);

% write out basin grid
fidw=fopen(BasinWinFile,'w');
fidp=fopen(BasinPreFile,'w');
fprintf(fidw,'%s\n',header);
fprintf(fidp,'%s\n',header);

for i=1:length(D.Basin.time)
    
   t=datestr(D.Basin.time(i),30);
   t([9 14 15])=[];
   header=sprintf(time_string,D.Basin.iLat(i),D.Basin.iLong(i),...
       D.Basin.DX(i),D.Basin.DY(i),D.Basin.SWLat(i),D.Basin.SWLon(i),t);
   fprintf('%s\n',header)
   
   out=D.Basin.Pre{i};
   fprintf(fidp,'%s\n',header);
   %fprintf('   Min,Max pre = %f %f\n',min(out(:)),max(out(:)))
   fprintf(fidp,value_string,out');
   fprintf(fidp,'\n');
   
   out=D.Basin.WinU{i};
   fprintf(fidw,'%s\n',header);
   %fprintf('   Min,Max u = %f %f\n',min(out(:)),max(out(:)))
   fprintf(fidw,value_string,out');
   fprintf(fidw,'\n');
   
   out=D.Basin.WinV{i};
   %fprintf('   Min,Max v = %f %f\n',min(out(:)),max(out(:)))
   fprintf(fidw,value_string,out');
   fprintf(fidw,'\n');
  
end

fclose(fidp);
fclose(fidw);

% write out region grid
if Region
    fidw=fopen(RegionWinFile,'w');
    fidp=fopen(RegionPreFile,'w');
    fprintf(fidw,'%s\n',header);
    fprintf(fidp,'%s\n',header);
    
    for i=1:length(D.Region.time)
        
        t=datestr(D.Region.time(i),30);
        t([9 14 15])=[];
        header=sprintf(time_string,D.Region.iLat(i),D.Region.iLong(i),...
            D.Region.DX(i),D.Region.DY(i),D.Region.SWLat(i),D.Region.SWLon(i),t);
        fprintf('%s\n',header)
        
        out=D.Region.Pre{i};
        fprintf(fidp,'%s\n',header);
        %fprintf('   Min,Max pre = %f %f\n',min(out(:)),max(out(:)))
        fprintf(fidp,value_string,out');
        fprintf(fidp,'\n');
        
        out=D.Region.WinU{i};
        fprintf(fidw,'%s\n',header);
        %fprintf('   Min,Max u = %f %f\n',min(out(:)),max(out(:)))
        fprintf(fidw,value_string,out');
        fprintf(fidw,'\n');
        
        out=D.Region.WinV{i};
        %fprintf('   Min,Max v = %f %f\n',min(out(:)),max(out(:)))
        fprintf(fidw,value_string,out');
        fprintf(fidw,'\n');
        
    end
    
    fclose(fidp);
    fclose(fidw);
    
end

