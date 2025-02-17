function WriteOwi(D,filename)
%WriteOwi - output OWI-formatted wind/pressure files from an OwiStruct 
%
%  original code for Basin and Region domains          B Blanton
%  code for Local domain                               R Luettich    2/9/2022
%  bug fix in 1st header line in Regional, Local files R Luettich 2/14/2022
%  eliminated minutes in 1st header line, all files    R Luettich 2/28/2022
%  eliminated extra digit to right of decimal pt in SWLon  R Luettich 1/29/2025
%  eliminated blank line between time snaps when # output values 
%                                   is divisible by 8  R Luettich 1/31/2025
%
% WriteOwi(Owi,filename)
% 
% INPUT:    D        - an OwiStruct .  See OwiStruct for details.
%           filename - prefix for file set.  
%                      If empty, WriteOwi writes to fort.2{21,22,23,24,17,18} as needed. 
%                      If a string, WriteOwi appends "win", "pre", "_Basin",  
%                      "_Region", "_Local" as needed.
%                      To specify the filename for each OWI file, pass in a 
%                      cell with 6 values.  E.g.:
%                      {'hur_Basin.pre','hur_Basin.win','hur_Region.pre',
%                       'hur_Region.win','hur_Local.pre','hur_Local.win'}
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

Local=true;
if ~isfield(D,'Local') || isempty(D.Local)
    fprintf('No local specified.  So no 217,218 files.\n')
    Local=false;
else
    if ~all(isfield(D.Local,flds))
        error('Input OwiStruct.Local is missing required fields.')
    end
end

% specify output file

if ~exist('filename')
   filename='fort';
   fprintf('No output file name specified, so using "fort".\n')
   % exit if fort.22{1,2,3,4,5,6} files already exist.  Will NOT overwrite...
   if  any([exist([filename '.221']) exist([filename '.222'])  ...
            exist([filename '.223']) exist([filename '.224'])  ...
            exist([filename '.217']) exist([filename '.218'])])
       str={'At least one fort.22{1,2,3,4} or fort.21{7,8} file exists and I cowardly refuse to overwrite.'};
       error(str{:});
   end
   BasinPreFile ='fort.221'; 
   BasinWinFile ='fort.222';
   RegionPreFile='fort.223'; 
   RegionWinFile='fort.224';
   LocalPreFile ='fort.217'; 
   LocalWinFile ='fort.218';
else
    if ischar(filename)
        BasinPreFile=[filename '_Basin.pre'];
        BasinWinFile=[filename '_Basin.win'];
        RegionPreFile=[filename '_Region.pre'];
        RegionWinFile=[filename '_Region.win'];
        LocalPreFile=[filename '_Local.pre'];
        LocalWinFile=[filename '_Local.win'];
    elseif iscell(filename)
        if length(filename)==6
            BasinPreFile=filename{1};
            BasinWinFile=filename{2};
            RegionPreFile=filename{3};
            RegionWinFile=filename{4};
            LocalPreFile=filename{5};
            LocalWinFile=filename{6};
        else
            error('Cant make sense of filename input.  Must be length==6 if a cell.')
        end
    end
end

%  set up the output formats

%time_string
%='iLat=%4diLong=%4dDX=%6.4fDY=%6.4fSWLat=%8.5fSWLon=%8.4fDT=%12s'; %original
time_string ='iLat=%4diLong=%4dDX=%6.4fDY=%6.4fSWLat=%8.5fSWLon=%8.3fDT=%12s';
value_string=' %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n';

% set Basin times
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

% compare times if Local
if Local
    lt0=D.Local.time(1);
    lt1=D.Local.time(end);
    if (bt0 ~= lt0) && (bt1~=lt1)
        error('Time span between Basin and Local are not equal.')
    end
    if length(D.Basin.time) ~= length(D.Local.time)
        error('Number of Time levels between Basin and Local is not equal.')
    end
end

t1=datestr(bt0,30);
t1([9 12 13 14 15])=[];
t2=datestr(bt1,30);
t2([9 12 13 14 15])=[];
header1=sprintf('Oceanweather WIN/PRE Format                            %10s     %10s',t1,t2);

% write out Basin files
fprintf('Writing Basin files:  ')
fidw=fopen(BasinWinFile,'w');
fidp=fopen(BasinPreFile,'w');
fprintf(fidw,'%s\n',header1);
fprintf(fidp,'%s\n',header1);

nextline_Basin=mod(D.Basin.iLat(1)*D.Basin.iLong(1),8);

for i=1:length(D.Basin.time)
    
   t=datestr(D.Basin.time(i),30);  % 30 stand for the format yyyymmddTHHMMSS
   t([9 14 15])=[];                % elminate T and SS
   header=sprintf(time_string,D.Basin.iLat(i),D.Basin.iLong(i),...
       D.Basin.DX(i),D.Basin.DY(i),D.Basin.SWLat(i),D.Basin.SWLon(i),t);
   %fprintf('%s\n',header);
   
   out=D.Basin.Pre{i};
   fprintf(fidp,'%s\n',header);   
   %fprintf('   Min,Max pre = %f %f\n',min(out(:)),max(out(:)))
   fprintf(fidp,value_string,out');
   if nextline_Basin~=0
       fprintf(fidp,'\n');
   end
   
   out=D.Basin.WinU{i};
   fprintf(fidw,'%s\n',header);
   %fprintf('   Min,Max u = %f %f\n',min(out(:)),max(out(:)))
   fprintf(fidw,value_string,out');
   if nextline_Basin~=0
       fprintf(fidw,'\n');
   end
   
   out=D.Basin.WinV{i};
   %fprintf('   Min,Max v = %f %f\n',min(out(:)),max(out(:)))
   fprintf(fidw,value_string,out');
   if nextline_Basin~=0   
       fprintf(fidw,'\n');
   end
  
end

fclose(fidp);
fclose(fidw);
fprintf('Done..\n')

% write out Region files
if Region
    fprintf('Writing Region files:  ')
    fidw=fopen(RegionWinFile,'w');
    fidp=fopen(RegionPreFile,'w');
    fprintf(fidw,'%s\n',header1);
    fprintf(fidp,'%s\n',header1);

    nextline_Region=mod(D.Region.iLat(1)*D.Region.iLong(1),8);    
    
    for i=1:length(D.Region.time)
        
        t=datestr(D.Region.time(i),30);
        t([9 14 15])=[];
        header=sprintf(time_string,D.Region.iLat(i),D.Region.iLong(i),...
            D.Region.DX(i),D.Region.DY(i),D.Region.SWLat(i),D.Region.SWLon(i),t);
        %fprintf('%s\n',header);
        
        out=D.Region.Pre{i};
        fprintf(fidp,'%s\n',header);
        %fprintf('   Min,Max pre = %f %f\n',min(out(:)),max(out(:)))
        fprintf(fidp,value_string,out');
        if nextline_Region~=0        
            fprintf(fidp,'\n');
        end
        
        out=D.Region.WinU{i};
        fprintf(fidw,'%s\n',header);
        %fprintf('   Min,Max u = %f %f\n',min(out(:)),max(out(:)))
        fprintf(fidw,value_string,out');
        if nextline_Region~=0         
            fprintf(fidw,'\n');
        end
        
        out=D.Region.WinV{i};
        %fprintf('   Min,Max v = %f %f\n',min(out(:)),max(out(:)))
        fprintf(fidw,value_string,out');
        if nextline_Region~=0         
            fprintf(fidw,'\n');
        end
        
    end
    
    fclose(fidp);
    fclose(fidw);
    fprintf('Done..\n')
end

% write out Local files
if Local
    fprintf('Writing Local files:  ')
    fidw=fopen(LocalWinFile,'w');
    fidp=fopen(LocalPreFile,'w');
    fprintf(fidw,'%s\n',header1);
    fprintf(fidp,'%s\n',header1);

    nextline_Local=mod(D.Local.iLat(1)*D.Local.iLong(1),8);  
    
    for i=1:length(D.Local.time)
        
        t=datestr(D.Local.time(i),30);
        t([9 14 15])=[];
        header=sprintf(time_string,D.Local.iLat(i),D.Local.iLong(i),...
            D.Local.DX(i),D.Local.DY(i),D.Local.SWLat(i),D.Local.SWLon(i),t);
        %fprintf('%s\n',header);
        
        out=D.Local.Pre{i};
        fprintf(fidp,'%s\n',header);
        %fprintf('   Min,Max pre = %f %f\n',min(out(:)),max(out(:)))
        fprintf(fidp,value_string,out');
        if nextline_Local~=0         
            fprintf(fidp,'\n');
        end
        
        out=D.Local.WinU{i};
        fprintf(fidw,'%s\n',header);
        %fprintf('   Min,Max u = %f %f\n',min(out(:)),max(out(:)))
        fprintf(fidw,value_string,out');
        if nextline_Local~=0         
            fprintf(fidw,'\n');
        end
        
        out=D.Local.WinV{i};
        %fprintf('   Min,Max v = %f %f\n',min(out(:)),max(out(:)))
        fprintf(fidw,value_string,out');
        if nextline_Local~=0 
            fprintf(fidw,'\n');
        end
        
    end
    
    fclose(fidp);
    fclose(fidw);
    fprintf('Done..\n')
end
