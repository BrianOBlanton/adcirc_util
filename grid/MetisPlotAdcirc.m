function [g,h]=MetisPlotAdcirc(PeStart,PeEnd)
%MetisPlotAdcirc
%Call as: [g,h]=MetisPlotAdcirc(PeStart,PeEnd)
%OR:      MetisPlotAdcirc(MetisStruct)
%
%See also: MetisLoadPeStruct, MetisShowPe

verbose=false;

DrawOnly=false;
if exist('PeStart','var')
   if isa(PeStart,'cell')   
      DrawOnly=true;
      MetisStruct=PeStart;
   end
end

cols=get(gca,'ColorOrder');
%eliminate black from color order
iding=all(cols'==0);
cols(iding,:)=[];

if ~DrawOnly
   % get sub-domain grids
   PEdirs=dir('PE*');
   if isempty(PEdirs)
      error('No metis domain decomposition directories found.')
   end

   if nargin==0
      PeStart=0;
      PeEnd=length(PEdirs)-1;
   elseif nargin==1
      PeEnd=PeStart;
   elseif nargin==2
      if PeStart<-1
         error('PeStart must be >=0')
      end
   end
   if PeStart>PeEnd
      error('PeStart must be <= PeEnd')
   end
   if PeEnd>=length(PEdirs)
      error('PeEnd must be < total number of subdomains.')
   end
   sprintf('Loading pes from PeStart=%d to PeEnd=%d',PeStart,PeEnd)

   MetisStruct=MetisLoadPeStruct(PeStart,PeEnd);
   
else
    
   for i=1:length(MetisStruct) 
      ii=rem(i,length(cols))+1;
      g=MetisStruct{i};
      h(i)=plotbnd(g,'Color',cols(ii,:),'LineWidth',1.5);
      %h(i)=drawelems(g,'Color',cols(ii,:),'LineWidth',.25);
      text(mean([min(g.x) max(g.x)]),mean([min(g.y) max(g.y)]),2,sprintf('%04d',g.Pe),...
         'HorizontalAlignment','center','VerticalAlignment','middle');
   end
   
end

%h(1)=plotbnd(gm,'Color','k','LineWidth',1);

drawnow
