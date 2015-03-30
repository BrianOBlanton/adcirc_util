function h=labcont(varargin)
%LABCONT label contours drawn with LCONTOUR
%
%  LABCONT labels contours on the current axes that have been
%  previously drawn with LCONTOUR.  
%  
%  LABCONT(p1,v1,p2,v2,...) passes the pv argument pairs to the 
%  text command for controlling the text properties. click on the contour 
%  to label with the left mouse-button. Continue clicking with the left 
%  mouse-button button until all the labels have been added, then press 
%  the right mouse-button in the current axes to stop LABCONT. 
%  
%  LABCONT('auto',ContourHandles,stride,p1,v1,p2,v2,...) puts text
%  labels on the contours referenced in ContourHandles, with a point stride
%  of stride. ContourHandles is a vector of graphics handles as returned by
%  LCONTOUR.  The auto-labeling is relatively crude in that LABCONT does
%  not try to optimize the lcoation of the text; use stride to coarsen the
%  text placement. 
%
%  If an output argument was supplied to LABCONT, the handles
%  to the text objects pointing to the labels is returned
%  to the MATLAB workspace.
%
% Written by : Brian O. Blanton, March 96
% added auto draft labeling, Nov 2013
%



%
% save the current value of the current figure's WindowButtonDownFcn,
% WindowButtonMotionFcn, and WindowButtonUpFcn
%
WindowButtonDownFcn=get(gcf,'WindowButtonDownFcn');
WindowButtonMotionFcn=get(gcf,'WindowButtonMotionFcn');
WindowButtonUpFcn=get(gcf,'WindowButtonUpFcn');
set(gcf,'WindowButtonDownFcn','');
set(gcf,'WindowButtonMotionFcn','');
set(gcf,'WindowButtonUpFcn','');


% if the first argin is 'auto', then auto label the axes
if nargin > 0 && strcmpi(varargin{1},'auto')
        varargin(1)=[];
        % second argument must be the contour handles
        hc=varargin{1};
        varargin(1)=[];
        hc(hc==0)=[];
        % third argument must be the point stride
        stride=varargin{1};
        if ~isa(stride,'double')
            disp({'Stride argument to LABCONT must be an integer.','LABCONT(''auto'',ContourHandles,stride,p1,v1,p2,v2,...)'})
            return;
        end
        varargin(1)=[];
        if stride<50, % warn user
            disp('A stride of < 50 may cause way too many text objects to be drawn.  Consider a larger number.')
        end
        
        for i=1:length(hc)
            t=get(hc(i),'UserData');
            cval=get(hc(i),'Color');

            xc=get(hc(i),'XData');
            yc=get(hc(i),'YData');
            xc=xc(1:stride:end);
            yc=yc(1:stride:end);
            t=repmat(int2str(t),[length(xc) 1]);
            ht{i}=text(xc,yc,t,'HorizontalAlignment','Center','Color',cval,varargin{:});
        end
        ht=cell2mat(ht(:));
else
    
    npts=0;
    fprintf('Click on contour with left mouse-button, right button to end.\n')

    while 1
       waitforbuttonpress;
       seltype=get(gcf,'SelectionType');
       if ~strcmp(seltype,'normal'),break,end

       target=gco;
       if ~strcmp(get(target,'Tag'),'contour')
          disp('Target NOT a contour from LCONTOUR or ISOPHASE')
       else
          npts=npts+1;
          val=get(target,'UserData');
          cval=get(target,'Color');
          pt=gcp;
          ht(npts)=text(pt(2),pt(4),num2str(val),...
                       'HorizontalAlignment','Center','Color',cval,'EdgeColor',cval,'BackgroundColor','w',varargin{:});
       end
    end

end

if nargout==1
   h=ht;
end

%
% return the saved values of the current figure's WindowButtonDownFcn,
% WindowButtonMotionFcn, and WindowButtonUpFcn to the current figure
%
set(gcf,'WindowButtonDownFcn',WindowButtonDownFcn);
set(gcf,'WindowButtonMotionFcn',WindowButtonMotionFcn);
set(gcf,'WindowButtonUpFcn',WindowButtonUpFcn);

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
