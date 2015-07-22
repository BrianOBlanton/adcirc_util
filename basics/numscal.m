function h=numscal(varargin)
%NUMSCAL label scalar (Q) values at x,y locations on current axes.
%
% Input:  x  - x-coordinate array
%         y  - y-coordinate array
%         Q  - scalar (1-D vector)
%
% NOTES: 1) NUMSCAL overlays existing plots regardless of
%           the state of hold
%        2) NUMSCAL determines if the current density of nodes
%           is too high to make the node numbers readable.
%           If node density is too high, NUMSCAL displays 
%           a warning box asking the user if NUMSCAL
%           should plot the numbers anyway.
%
% Call as: >> h=numscal(x,y,Q,p1,v1,...);

%
% Written by : Brian O. Blanton
%

% DEFINE ERROR STRINGS
err1=['Not enough input arguments; need atleast x,y,q or fgs,q'];
err2=['Lengths of x,y must be the same'];
err3=['Size of Q must be the same as x,y'];
warn1=str2mat('The current axes is too dense for the  ',...
       'scalar values to be readable.  CONTINUE?'); 
   
% check arguments
if nargin==0
   disp('h=numscal(x,y,Q,p1,v1,...);');
   return
end
if nargin < 2 
   error(err1);
end

if isstruct(varargin{1})
   fgs=varargin{1};
   varargin(1)=[];
   x=fgs.x;
   y=fgs.y;
   Q=varargin{1};
   varargin(1)=[];
   if ischar(Q)
       if strcmp(Q,'z')
           Q=fgs.z;
       else
           error('If arg 2|3 is a string, it must be ''z''')
       end
   end 

else
   x=varargin{1};
   y=varargin{2};
   Q=varargin{3};
   varargin([1 2 3])=[];
end

if length(Q)>length(x)
    % must be element-based
    if ~isfield(fgs,'xecen')
        error({'If length Q == length(e(:,1), then xecen, yecen must be fields in fem_grid_struct.'})
    end
    x=fgs.xecen;
    y=fgs.yecen;
end

% Default propertyname values
Digits=6;

% Strip off propertyname/value pairs in varargin not related to
% "line" object properties.
k=1;
while k<length(varargin),
  switch lower(varargin{k}),
    case 'digits',
      Digits=varargin{k+1};
      varargin([k k+1])=[];
    otherwise
      k=k+2;
  end;
end;

if length(varargin)<2
   varargin={};
end

% check input arguments length
if length(x)~=length(y) 
   error(err2);
end
if length(x)~=length(Q) 
   error(err3);
end

% compute offset as 2% of each axis range
X=get(gca,'Xlim');
Y=get(gca,'YLim');
%off=(X(2)-X(1))/50;
off=0;

% get indices of nodes within viewing window defined by X,Y
filt=find(x>=X(1)&x<=X(2)&y>=Y(1)&y<=Y(2));

% determine if viewing window is zoomed-in
% enough for node numbers to be meaningful
oldunits=get(gca,'Units'); 
set(gca,'Units','Points'); 
rect=get(gca,'Position');        % get viewing window width in 
set(gca,'Units',oldunits);       % point units (1/72 inches)
xr=rect(3)-rect(1);
yr=rect(4)-rect(2);
xden=xr/sqrt(length(filt));
yden=yr/sqrt(length(filt));
%den=sqrt(xden*xden+yden*yden);
%if den < 5*ps
%   click=questdlg(warn1,'yes','no');
%   if strcmp(click,'no'),return,end
%end


% label only those nodes that lie within viewing window.
%   h1=text(x(filt)+off,y(filt)+off,int2str(Q(filt)),...
fmtstr=sprintf('%%.%dg',Digits);
   h1=text(x(filt)+off,y(filt)+off,num2str(Q(filt),fmtstr),...
       'HorizontalAlignment','center',...
       'VerticalAlignment','middle',...
       'Color','k',...
       'EdgeColor','k',...
       'BackgroundColor','w',...
       'Tag','Node Scalar Value',...
       'FontSize',16,varargin{:});
        
if nargout>0,h=h1;end


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



