function retval=plotbndadc(fem_grid_struct,varargin)
%PLOTBNDADC plot color-coded boundary of FEM mesh
%   PLOTBNDADC(fem_grid_struct) draws a color-coded boundary of
%   an ADCIRC FEM domain.
%
%
%   INPUT : fem_grid_struct (from LOADGRID, see FEM_GRID_STRUCT)
%           standard LINE parameter/value pairs can be passed in.
%
%  OUTPUT : hboun - handle to boundary object drawn
%
%    CALL : hboun=plotbnd(fem_grid_struct,p1,v1,...)
%
% Written by : Brian O. Blanton
% Summer 2007
%

% added varargin, 28 Jun 02, BOB

%if nargin ~=1
%   error('    Incorrect number of input arguments to PLOTBND');
%end

if ~is_valid_struct(fem_grid_struct)
    error('    Argument to PLOTBNDADC must be a valid fem_grid_struct.')
end

% Extract grid fields from fem_grid_struct
%
bnd=fem_grid_struct.bnd;
x=fem_grid_struct.x;
y=fem_grid_struct.y;

% this is the complete boundary
ns=bnd(:,1);
ne=bnd(:,2);
X=[x(ns) x(ne) NaN*ones(size(ns))]';
Y=[y(ns) y(ne) NaN*ones(size(ns))]';
X=X(:);
Y=Y(:);
hboun=line(X,Y,'Tag','boundary','Color','k',varargin{:});

% legstr={'Land'};
% hleg=hboun;

legstr={};
hleg=[];

cols={
    'g' 'none' '*'   % type 0
    'm' 'none' '.'   % type 1
    '-' 'none' '-'   % type 2
    '-' 'none' '-'   % type 3
    '-' 'none' '-'   % type 4
    '-' 'none' '-'   % type 5
    '-' 'none' '-'   % type 6 (n/a)
    '-' 'none' '-'   % type 7 (n/a)
    '-' 'none' '-'   % type 8 (n/a)
    '-' 'none' '-'   % type 9 (n/a)
    'b' 'none' '.'   % type 10
    'c' 'none' '.'   % type 11
    'r' 'none' '.'   % type 12
    'c' 'none' '.'   % type 13
    '-' 'none' '-'   % type 14 (n/a)
    '-' 'none' '-'   % type 15 (n/a)
    '-' 'none' '-'   % type 16 (n/a)
    '-' 'none' '-'   % type 17 (n/a)
    '-' 'none' '-'   % type 18 (n/a)
    '-' 'none' '-'   % type 19 (n/a)
    'c' 'none' '.'   % type 20
    'b' 'none' '.'   % type 21
    'r' 'none' '.'   % type 22
    'c' 'none' '.'   % type 23
    'm' 'none' 's'   % type 24
    'k' 'none' 's'   % type 64 VEW
    };

legs={'Land (0, strong no normal, free tangential)'            % ibtype==0
      'Island (1, strong no normal, free tangential)'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      'Land (10, strong no normal, no tangential)'
      'Island (11, strong no normal, no tangential)'
      'River (12, strong normal, no tangential)'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      '-'
      'Land (20, weak no normal, free tangential)'
      'Island (21, weak no normal, free tangential)'
      'River (22, weak normal, free tangential)'
      '-'
      'Weir (24)'
      'VEW (64)'
      };

nriver=0;
if isfield(fem_grid_struct,'ibtype')
    
    % Boundary types (not open) in grid
    Ibtypes=unique(fem_grid_struct.ibtype);
    
    for i=1:length(Ibtypes)
        
        % +1 bec ibtype starts at 0 and the cols and 
        % legs cells are indexed by ibtype
        ibtype=Ibtypes(i)+1;  
        idx=find(fem_grid_struct.ibtype==Ibtypes(i));

        fprintf(' Drawing ibtype=%d, %d segments\n',Ibtypes(i),length(idx))

        switch Ibtypes(i)
            
            case 0  % land (no normal, free tangential)
                
                for k=1:length(idx)
                    j=fem_grid_struct.ln{idx(k)};
                    x=fem_grid_struct.x(j);
                    y=fem_grid_struct.y(j);
                    h=line(x,y,'Color',cols{ibtype,1},'LineStyle','-','Marker','.'); 
                    %line(x(1:10:end),y(1:10:end),'Color',cols{ibtype,1},'LineStyle','non','Marker',cols{ibtype,3},'MarkerSize',8);
                    % temp=repmat(ibtype-1,[length(x) 1]);
                    % text(x(1:10:end),y(1:10:end),string(temp(1:10:end)),'Color',cols{ibtype,1},'FontSize',10);

                    % text(x(floor(length(x)/2)),y(floor(length(x)/2)),int2str(idx(k)),'EdgeColor',cols{ibtype,1},'BackgroundColor','w')
                end
                hleg=[hleg;h];
                legstr={legstr{:},legs{ibtype}};
                
            case 1  % island
                
                for j=1:length(idx)
                    k=fem_grid_struct.ln{idx(j)};
                    x=fem_grid_struct.x(k);
                    y=fem_grid_struct.y(k);
                    h=line(x,y,'Color',cols{jj,1},'LineStyle',cols{jj,2},'Marker',cols{jj,3},'MarkerSize',12);
                end
                hleg=[hleg;h];
                legstr={legstr{:},legs{jj}};
                
            case 10  % Land (10, strong no normal, no tangential)
                
                for j=1:length(idx)
                    k=fem_grid_struct.ln{idx(j)};
                    x=fem_grid_struct.x(k);
                    y=fem_grid_struct.y(k);
                    h=line(x,y,'Color',cols{jj,1},'LineStyle',cols{jj,2},'Marker',cols{jj,3},'MarkerSize',12);
                end
                hleg=[hleg;h];
                legstr={legstr{:},legs{jj}};
                
            case 12  % River (12, strong normal, no tangential)
                
                for j=1:length(idx)
                    k=fem_grid_struct.ln{idx(j)};
                    x=fem_grid_struct.x(k);
                    y=fem_grid_struct.y(k);
                    h=line(x,y,'Color',cols{jj,1},'LineStyle',cols{jj,2},'Marker',cols{jj,3},'MarkerSize',12);
                end
                hleg=[hleg;h];
                legstr={legstr{:},legs{jj}};
                
                
            case 20  % Land (20, weak no normal, free tangential)
                
                for j=1:length(idx)
                    k=fem_grid_struct.ln{idx(j)};
                    x=fem_grid_struct.x(k);
                    y=fem_grid_struct.y(k);
                    h=line(x,y,'Color',cols{ibtype,1},'LineStyle','-','Marker','.');                    
                    
                    % temp=repmat(ibtype-1,[length(x) 1]);
                    % text(x(1:10:end),y(1:10:end),string(temp(1:10:end)),'Color',cols{ibtype,1},'FontSize',10);

                    %line(x(1:10:end),y(1:10:end),'Color',cols{ibtype,1},'LineStyle','none','Marker',cols{ibtype,3},'MarkerSize',12);
                    text(x(floor(length(x)/2)),y(floor(length(x)/2)),int2str(idx(j)),...
                        'Color',cols{ibtype,1},'EdgeColor',cols{ibtype,1},'BackgroundColor','w')
                end
                hleg=[hleg;h];
                legstr={legstr{:},legs{ibtype}};
                
            case 21  % Island (21, weak no normal, free tangential)
                
                for j=1:length(idx)
                    k=fem_grid_struct.ln{idx(j)};
                    x=fem_grid_struct.x(k);
                    y=fem_grid_struct.y(k);
                    h=line(x,y,'Color',cols{ibtype,1},'LineStyle','-','Marker','.');
%                    line(x(1:10:end),y(1:10:end),'Color',cols{ibtype,1},'LineStyle','none','Marker',cols{ibtype,3},'MarkerSize',12);
                    % temp=repmat(ibtype-1,[length(x) 1]);
                    % text(x(1:10:end),y(1:10:end),string(temp(1:10:end)),'Color',cols{ibtype,1},'FontSize',10);
                    %text(mean(x),mean(y),int2str(idx(j)),'EdgeColor',cols{ibtype,1},'BackgroundColor','w')
                    text(x(floor(length(x)/2)),y(floor(length(x)/2)),int2str(idx(j)),...
                        'Color',cols{ibtype,1},'EdgeColor',cols{ibtype,1},'BackgroundColor','w')
                end
                hleg=[hleg;h];
                legstr={legstr{:},legs{ibtype}};
                
            case 22  % River (22, weak normal, free tangential)
                              
                for j=1:length(idx)
                    nriver=nriver+1;
                    k=fem_grid_struct.ln{idx(j)};
                    x=fem_grid_struct.x(k);
                    y=fem_grid_struct.y(k);
                    h=line(x,y,'Color',cols{jj,1},'LineStyle',cols{jj,2},'Marker',cols{jj,3},'MarkerSize',12);
                    text(mean(x),mean(y),int2str(nriver),'EdgeColor',cols{jj,1},'BackgroundColor','w')
                end
                hleg=[hleg;h];
                legstr={legstr{:},legs{jj}};
                
            case 24
                i24=find(fem_grid_struct.ibtype==24);
                for ii=1:length(i24)
                    j=fem_grid_struct.ln{i24(ii)};
                    x=fem_grid_struct.x(j);
                    y=fem_grid_struct.y(j);
                    %line(x,y,'Color','r',LineWidth',2)
                    %h=line(x,y,'Color',cols{jj,1},'LineStyle',cols{jj,2},'Marker',cols{jj,3},'MarkerSize',12);
                    h=line(x,y,'Color',cols{jj,1},'LineStyle',cols{jj,2},'Marker',cols{jj,3},'MarkerSize',6);

                    %   j=floor(length(x)/2);
                    %   text(x(j),y(j),int2str(i),'Horiz','center','Vert','middle')
                    
%                     j=fem_grid_struct.ln{i24(i)};
%                     x=fem_grid_struct.x(j);
%                     y=fem_grid_struct.y(j);
%                     h=line(x,y,'Color',cols{jj,1},'LineStyle',cols{jj,2},'Marker',cols{jj,3},'MarkerSize',12);
                end
                hleg=[hleg;h(1)];
                legstr={legstr{:},legs{jj}};
            otherwise
                
        end
        
    end
end

% color bnd types
if isfield(fem_grid_struct,'nopenboundaries')
    if fem_grid_struct.nopenboundaries>0
        fprintf(' Drawing elevation boundaries, %d segments\n',fem_grid_struct.nopenboundaries)
        for i=1:fem_grid_struct.nopenboundaries
            iob=fem_grid_struct.ob{i};
            x=fem_grid_struct.x(iob);
            y=fem_grid_struct.y(iob);
            hob=line(x,y,'Color','r','LineStyle','-','Marker','.');
            %line(x(1:5:end),y(1:5:end),'Color','r','LineStyle','none','Marker','x','MarkerSize',6);
            text(x(floor(length(x)/2)),y(floor(length(x)/2)),int2str(i),...
                'Color','r','EdgeColor','r','BackgroundColor','w')
        end
        legstr={legstr{:},'Elevation'};
        hleg=[hleg;hob(1)];
    end
end

if ~isempty(legstr)
    legend(hleg,legstr,'Location','NorthWest')
end

if nargout==1,retval=hboun;end
%
%LabSig  Brian O. Blanton
%        Department of Marine Sciences
%        12-7 Venable Hall
%        CB# 3300
%        University of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        brian_blanton@unc.edu
%
%        Summer 1997
%
