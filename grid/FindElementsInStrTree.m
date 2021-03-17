function j=FindElementsInStrTree(fgs,points_x,points_y,TOL)
% Call as:  j=FindElementsInStrTree(fgs,points_x,points_y);

[m,n]=size(points_x);
j=NaN*ones(m,n);
if ~exist('TOL')
    TOL=1e-2;
end

strtree=fgs.strtree;

%%
%tic;
for i=1:m

%parfor i=1:m
    
    %fprintf('%04d/%05d ... \n',i,m)
    temppx=points_x(i,:);
    temppy=points_y(i,:);
    jtemp=j(i,:);
    
    for ii=1:n
        px=temppx(ii);
        py=temppy(ii);
        p=com.vividsolutions.jts.geom.Coordinate(px,py);
        e=com.vividsolutions.jts.geom.Envelope(p);
        l=strtree.query(e);

        if l.size==0
            %fprintf('no potential elements found for %f, %f\n',px,py);
        else
            
%             v=NaN*ones(l.size,1);
%             for k=0:l.size-1
%                 v(k+1)=l.get(k);
%             end

            v=ArrayListToVector(l);
            phi=basis2d(fgs,[px*ones(size(v)) py*ones(size(v))],v)';
            idx=find(all(phi<=1+TOL & phi>=0-TOL));

            if ~isempty(idx)
%                j(i,ii)=v(idx(1));
                jtemp(ii)=v(idx(1));
            end

        end
    end
    j(i,:)=jtemp;
end

%t=toc;
%disp(sprintf('Containing elements for %d point(s) found in %.1f secs',m*n,t));

