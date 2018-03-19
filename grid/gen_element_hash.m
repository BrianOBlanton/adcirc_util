function fgsout=gen_element_hash(fgs,n,bbox)
%fgsout=gen_element_hash(fgs,n,box)
%
% box [xmin xmax ymin ymax], as in axis

debug=1;

if n<2
   error('N must be > 2 in GEN_ELEMENT_HASH');
end
ehash.n=n;

% get coords of element centroids
gtemp=attach_elem_centroids(fgs);

if exist('bbox')

   idx=find(gtemp.xecen>=bbox(1) & gtemp.xecen<=bbox(2) & ...
            gtemp.yecen>=bbox(3) & gtemp.yecen<=bbox(4)); 
   dx=fgs.A(idx,:);dx=dx(:);
   dy=fgs.B(idx,:);dy=dy(:);
   dd=sqrt(dx.*dx+dy.*dy);
   maxb=max(dd);
   minb=min(dd);
   
   ehash.x=linspace(bbox(1)-maxb,bbox(2)+maxb,ehash.n);
   ehash.y=linspace(bbox(3)-maxb,bbox(4)+maxb,ehash.n);
   ehash.dx2=ehash.x(2)-ehash.x(1);
   ehash.dy2=ehash.y(2)-ehash.y(1);

else
   
   dx=fgs.A(:);
   dy=fgs.B(:);
   dd=sqrt(dx.*dx+dy.*dy);
   maxb=max(dd);
   minb=min(dd);
 
   ehash.x=linspace(min(fgs.x)-maxb,max(fgs.x)+maxb,ehash.n);
   ehash.y=linspace(min(fgs.y)-maxb,max(fgs.y)+maxb,ehash.n);
   ehash.dx2=ehash.x(2)-ehash.x(1);
   ehash.dy2=ehash.y(2)-ehash.y(1);
   
end

dfuz=maxb/5;
%dfuz=mean(dd)*2;

if debug
   figure
   plotbnd(fgs)
   axis('equal')
   if exist('bbox')
      axis(bbox([1 2 3 4])+[-1 1 -1 1]*maxb)
   end
   line(repmat(ehash.x,[length(ehash.y) 1]),repmat(ehash.y',[1 length(ehash.x)]),'Color','k')
   line(repmat(ehash.x,[length(ehash.y) 1])',repmat(ehash.y',[1 length(ehash.x)])','Color','k')
   drawnow
end

% elements
jtemp=ceil((gtemp.xecen-ehash.x(1))/ehash.dx2);
itemp=ceil((gtemp.yecen-ehash.y(1))/ehash.dy2);

k=0;
Color={'r','g','b','m'};
ncol=length(Color);

for i=1:ehash.n-1
for j=1:ehash.n-1

   idx=find(itemp==i);
   jdx=find(jtemp==j);
   test1=intersect(idx,jdx);

   % add layer of elements surrounding each square.
   jdx=find(gtemp.xecen>ehash.x(j)-dfuz & gtemp.xecen<ehash.x(j+1)+dfuz);
   idx=find(gtemp.yecen>ehash.y(i)-dfuz & gtemp.yecen<ehash.y(i+1)+dfuz);
   test2=intersect(idx,jdx);
   test=union(test1,test2);
   %test=test1;
   if ~isempty(test)
      k=k+1;
      ehash.e{i,j}=test;
      if debug
%if k==10,keyboard,end
         line(gtemp.xecen(test),gtemp.yecen(test),'LineStyle','none','Color',Color{rem(k,ncol)+1},'Marker','.','MarkerSize',3);
         drawnow
      end
   end
end
end

fgsout=fgs;
fgsout.ehash=ehash;
com=sprintf('save %s.ehash.mat  ehash',fgs.name);
eval(com);
