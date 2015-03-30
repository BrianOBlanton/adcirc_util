function h=MetisShowPe(MetisStruct,PeStart,PeEnd)
% MetisShowPe(MetisStruct,PeStart,PeEnd)

cols=get(gca,'ColorOrder');
%eliminate black from color order
iding=all(cols'==0);
cols(iding,:)=[];

if nargin==0
    disp('h=MetisShowPe(MetisStruct,PeStart,PeEnd);');
    return
else
    if ~isa(MetisStruct,'cell')
        error('First argument to MetisShowPe must be a MetisStruct')
    end    
end

if nargin==1
    PeStart=0;
    PeEnd=length(MetisStruct)-1;
    PeList=(PeStart:PeEnd);
elseif nargin==2
    PeList=PeStart;
    PeStart=PeList(1);
    PeEnd=PeList(end);
elseif nargin==3
    if PeStart<-1
        error('PeStart must be >=0')
    end
    PeList=(PeStart:PeEnd);
end

if PeStart>PeEnd
   error('PeStart must be <= PeEnd')
end
if PeEnd>=length(MetisStruct)
   error('PeEnd must be < total number of subdomains.')
end

%[PeStart PeEnd]

disp(['Drawing pes ',int2str(PeList)])

j=0;
for i=PeList
    j=j+1;
    g=MetisStruct{i+1};
    ii=rem(j,length(cols))+1;
    h(j)=plotbnd(g,'Color',cols(ii,:));
    %h(i+1)=drawelems(g,'Color',cols(ii,:),'LineWidth',.25);
    text(mean([min(g.x) max(g.x)]),mean([min(g.y) max(g.y)]),2,sprintf('%04d',g.Pe),...
         'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','Color',cols(ii,:));
end


