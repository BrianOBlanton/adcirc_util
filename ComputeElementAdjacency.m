function icee=ComputeElementAdjacency(fgs)

if ~isfield(fgs,'xecen')
    disp('Attaching element centroids ...')
    fgs=attach_elem_centroids(fgs);
end
if ~isfield(fgs,'strtree')
    error('Need to add STRtree field to fgs.  Call ComputeStrTree first.')
end

%i=[fgs.e(:,1);fgs.e(:,2);fgs.e(:,3)];
%j=[fgs.e(:,2);fgs.e(:,3);fgs.e(:,1)];
ii=[fgs.e(:,1) fgs.e(:,2) fgs.e(:,3)];
jj=[fgs.e(:,2) fgs.e(:,3) fgs.e(:,1)];

%% Compute test points for evaluating element tests.

%% 1) compute line segment mid-points
midp.x=[fgs.x(ii)+(fgs.x(jj)-fgs.x(ii))/2];
midp.y=[fgs.y(ii)+(fgs.y(jj)-fgs.y(ii))/2];
%line(midp.x(:) ,midp.y(:),'Marker','+','Color','g','LineStyle','none')

%% 2) Compute boundary line segment slopes and replace Inf's by 1e10
diff= abs(fgs.x(jj)-fgs.x(ii));
numerator=fgs.y(jj)-fgs.y(ii);
denominator=fgs.x(jj)-fgs.x(ii);
idx=find(diff>eps);
slope=1e10*ones(size(diff));
slope(idx)=numerator(idx)./denominator(idx);

%% 3) Compute perpendicular slopes and replace Inf's by 1e10
slopeperp=-1./slope;
idx=find(slopeperp>1e10);
slopeperp(idx)=1e10*ones(size(idx));

%% 4) Compute unit direction vectors
mag=sqrt((fgs.x(jj)-fgs.x(ii)).*(fgs.x(jj)-fgs.x(ii))+(fgs.y(jj)-fgs.y(ii)).*(fgs.y(jj)-fgs.y(ii)));
UDV.x=[fgs.x(jj)-fgs.x(ii)]./mag;
UDV.y=[fgs.y(jj)-fgs.y(ii)]./mag;

%% 5) rotate vectors; theta=-pi/2, so cos terms==0, sin terms==-1
PUDV.x= UDV.y;
PUDV.y=-UDV.x;

%% 6) test points are along normal direction vectors eminating
%    from the boundary line-segment mid-points.
dxy=sqrt(2*fgs.ar)/10;
dxy=repmat(dxy,[1 3]);
%dxy=dxy(:);
test.x=midp.x+PUDV.x.*dxy;
test.y=midp.y+PUDV.y.*dxy;
%line(test.x,test.y,'Marker','o','Color','g','LineStyle','none')

disp(sprintf('Finding elements for %d test points ...',length(test.x)*3))
tic;j=FindElementsInStrTree(fgs,test.x,test.y);toc

icee=NaN*ones(fgs.ne,3);
tic;
for i=1:fgs.ne
   for j=1:3
   p=com.vividsolutions.jts.geom.Coordinate(test.x(i,j),test.y(i,j));   
   e=com.vividsolutions.jts.geom.Envelope(p);
   l=fgs.strtree.query(e);
   for k=0:l.size-1
      teste=l.get(k);
      phi=basis2d(fgs,[test.x(i,j) test.y(i,j)],teste);
      if all(phi<=1 & phi>=0),icee(i,j)=teste;break;end
   end
   end
end
t=toc;
disp(sprintf('Containing elements for %d points found in %.1f secs',fgs.ne*3,t));
