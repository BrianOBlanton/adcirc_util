function pes=MetisLoadPeStruct(PeStart,PeEnd)
%MetisLoadPeStruct
%Call as: pes=MetisLoadPeStruct(PeStart,PeEnd)
%
%See also: MetisPlotAdcirc, MetisShowPe

verbose=false;

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
fprintf('Loading pes from PeStart=%d to PeEnd=%d.\n',PeStart,PeEnd)

j=0; 
pes=cell(PeEnd-PeStart+1,1);

for i=PeStart:PeEnd

     j=j+1;
     fort14=[PEdirs(i+1).name '/fort.14'];
     if verbose
         fprintf('Loading %d(%d),%s ... \n',i,j,fort14)
     end
     if ~exist(fort14,'file')
        error(['fort.14 file for ' fort14 ' not found.'])
     end
     
     pes{j}=grd_to_opnml(fort14,verbose);
     pes{j}.Pe=i;
     
%         ii=rem(i,length(cols))+1;         
%         h(j)=plotbnd(g(j),'Color',cols(ii,:));
%         h(i+1)=drawelems(g(i),'Color',cols(ii,:),'LineWidth',.25);
%         text(mean([min(g(j).x) max(g(j).x)]),mean([min(g(j).y) max(g(j).y)]),2,sprintf('%04d',i),...
%            'HorizontalAlignment','center','VerticalAlignment','middle');

    % see if centroid is within bounding box
    x=mean(pes{j}.x);
    y=mean(pes{j}.y);
    xll=min(pes{j}.x);
    yll=min(pes{j}.y);
    xur=max(pes{j}.x);
    yur=max(pes{j}.y);
    inp=inpolygon(x,y,[xll xur xur xll xll],[yll yll yur yur yll]);
    pes{j}.split=~inp;
    if ~inp
        fprintf('%d pe might be split.\n',j);
    end

end
