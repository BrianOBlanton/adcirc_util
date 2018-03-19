function gen_sms_poly(fem_grid_struct)
%GEN_SMS_POLY generate an SMS poly and xyz file from a fem_grid_struct
% GEN_SMS_POLY(fem_grid_struct) generate an SMS poly and xyz 
% file from an existing fem_grid_struct.
%
% Written by : Brian O. Blanton and Alfredo Aretxabaleta
% Summer 2000


% Copy fem_grid_struct
gridstruct=fem_grid_struct;

if nargin ~=1
   error('    Incorrect number of input arguments to GEN_SMS_POLY');
end

% Invalid fem_grid_struct OK in this case. Make sure 
% .e,.x.,y,.z are there.
if ~isfield(fem_grid_struct,'e')
   disp('    Element list field not part of fem_grid_struct');return
elseif ~isfield(fem_grid_struct,'x')
   disp('    X-node list field not part of fem_grid_struct');return
elseif ~isfield(fem_grid_struct,'y')
   disp('    Y-node list field not part of fem_grid_struct');return
elseif ~isfield(fem_grid_struct,'z')
   disp('    Bathymetry list field not part of fem_grid_struct');return
end

if ~isfield(fem_grid_struct,'bnd')
   gridstruct.bnd=detbndy(fem_grid_struct.e);
   disp('Added .bnd to fem_grid_struct')
end
if ~isfield(fem_grid_struct,'ar')
   gridstruct=el_areas(gridstruct);
   disp('Added areas to fem_grid_struct')
end

   gridstruct=belint(gridstruct);



% Force orientation to be LHS
newbndlist=orient_bnd(gridstruct);

codes=NaN*ones(size(gridstruct.bnd(:,1)));

% Generate bel list for the exterior boundary
   % Starting node on bnd
   strnn=min(newbndlist(:));
   endnn=strnn;

   % Connect the Dots (Starting NN to Ending NN CCW!!)
   % Since the boundary is connected in CCW order, the list of boundary
   % segments must be oriented with the interior of the domain to the
   % left-hand side of the segment.  This is what 'HashBnd' does.
   %
   bnd=newbndlist;
   killbnd=newbndlist;
   connected=0;
   nextnn=strnn;
   disp('Connecting Start to End NN, CCW ...');
   code_lst=[];
   while ~connected
      connr=bnd(find(bnd(:,1)==nextnn),2);
      % Insert new code into bellist
      idx=find(bnd(:,1)==nextnn);
      code_lst=[code_lst
                idx];
      if connr==endnn,connected=1;,end
      nextnn=connr;
      idelete=find(killbnd(:,1)==nextnn);
      killbnd(idelete,:)=[];
   end

   % Set code list for land
   codes(1:length(code_lst))=1;
   bnd=bnd(code_lst,:);

   % The remainder of killbnd is all islands
   islandcount=0;
   if ~isempty(killbnd)
      finished=0;
      islbnd=killbnd;
      while ~finished
         islandcount=islandcount+1;
         disp(['Connecting Island ' int2str(islandcount) '...']);
	 % pick a node number
	 strnn=min(killbnd(:));
	 
	 endnn=strnn;
	 connected=0;
	 nextnn=strnn;
	 code_lst=[];
	 while ~connected
	    connr=islbnd(find(islbnd(:,1)==nextnn),2);
	    % Insert new code into bellist
	    idx=find(islbnd(:,1)==nextnn);
	    code_lst=[code_lst
                      idx];
	    if connr==endnn,connected=1;,end
	    nextnn=connr;
	    idelete=find(killbnd(:,1)==nextnn);
	    killbnd(idelete,:)=[];
	 end
	 inan=find(isnan(codes));
	 inan=inan(1);
	 codes(inan:inan+length(code_lst)-1)=islandcount+1;
	 bnd(inan:inan+length(code_lst)-1,:)=islbnd(code_lst,:);
	 if isempty(killbnd),finished=1;,end
      end
   end
   
blist=bnd(:,1);
nseg=length(unique(codes));
   
% Output SMS poly file
fid=fopen([gridstruct.name '.poly'],'w');
for i=1:nseg
   iseg=find(codes==i);
   
   iseglist=blist(iseg);
   
   out=[gridstruct.x(iseglist); gridstruct.y(iseglist); gridstruct.z(iseglist)];
   x=gridstruct.x(iseglist);
   y=gridstruct.y(iseglist); 
   z=gridstruct.z(iseglist);
   
   fprintf(fid,'POLY %d\n',length(iseg));
   fprintf(fid,'%f %f %f\n',[x y z]');
   
end

fclose(fid);

% Output SMS xyz file
fid=fopen([gridstruct.name '.xyz'],'w');
fprintf(fid,'XYZ\n');
fprintf(fid,'%f %f %f\n',[gridstruct.x gridstruct.y gridstruct.z]');



fclose(fid);


   
return



function bnd=orient_bnd(gridstruct)
   
   e=gridstruct.e;
   x=gridstruct.x;
   y=gridstruct.y;
   AR=gridstruct.ar;
   A=gridstruct.A;
   B=gridstruct.B;
   T=gridstruct.T;

   
   % Determine Boundary list
   %
   % Form (i,j) connection list from .ele element list
   %
   i=[e(:,1);e(:,2);e(:,3)];
   j=[e(:,2);e(:,3);e(:,1)];

   % Form the sparse adjacency matrix and add transpose.
   %
   n = max(max(i),max(j));
   ICM = sparse(i,j,-1,n,n);
   ICM = ICM + ICM';

   % Consider only the upper part of ICM, since ICM is symmetric
   % 
   ICM=ICM.*triu(ICM);

   % The boundary segments are ICM's with value == 1
   %
   ICM=ICM==1;

   % Extract the row,col from new ICM for the boundary list.
   %
   [ib,jb,s]=find(ICM);
   ib=ib(:);jb=jb(:);

   % Sort Col#1 of bnd and permute Col#2
   %
   [ib,iperm]=sort(ib);
   jb=jb(iperm);
   bnd=[ib(:) jb(:)];

   disp('Re-ordering boundary list (this may take a while) ...');
   
   % The boundary list generated by the above sparse-matrix method
   % does not ensure that the "left-hand toward the grid interior" 
   % convention is maintained.  We'll force that part now.  Unfortunately,
   % this could take a while.  However, once this is done for a mesh, it 
   % need not be done again.  An existing .bel file can be loaded; its boundary
   % list will necessarily be ordered correctly. 
   %

   % Compute test points as follows:
   %
   %      Assume that the segments are oriented such that
   %      the interior of the FEM domain is to the left and then:
   %         1) compute boundary line-segment mid-points
   %         2) compute slopes of boundary line-segments
   %         3) compute normal slopes of boundary line-segments
   %         4) compute unit direction vectors of the boundary segment
   %            from the first node to the second.
   %         5) rotate the unit vectors pi/2 CCW, which SHOULD point 
   %            toward the interior, atleast locally.
   %         6) test points lie along normal segment direction vectors
   %            from mid-point.

   % 1) compute boundary line segment mid-points
        midp=[x(ib)+(x(jb)-x(ib))/2 y(ib)+(y(jb)-y(ib))/2];

   % 2) Compute boundary line segment slopes and replace Inf's by 1e10
        diff= abs(x(jb)-x(ib));
	numerator=y(jb)-y(ib);
	denominator=x(jb)-x(ib);
	idx=find(diff>eps);
	slope=1e10*ones(size(diff));
	slope(idx)=numerator(idx)./denominator(idx);

   % 3) Compute perpendicular slopes and replace Inf's by 1e10
	slopeperp=-1./slope;
	idx=find(slopeperp>1e10);
	slopeperp(idx)=1e10*ones(size(idx));

   % 4) Compute unit direction vectors 
	mag=sqrt((x(jb)-x(ib)).*(x(jb)-x(ib))+(y(jb)-y(ib)).*(y(jb)-y(ib)));
	UDVx=[x(jb)-x(ib)]./mag;
	UDVy=[y(jb)-y(ib)]./mag;

   % 5) rotate vectors; theta=pi/2, so cos terms==0, sin terms==1
	PUDVx=-UDVy;
	PUDVy=UDVx;

   % 6) test points are along normal direction vectors eminating
   %    from the boundary line-segment mid-points.
	dxy=sqrt(min(AR))/10;
	testx=midp(:,1)+PUDVx*dxy;
	testy=midp(:,2)+PUDVy*dxy;

   disp('Ignore "Divide by zero" warnings!!');
   
   % Now, locate elements for each test point.  If a point is in the 
   % domain, an element will be found for it, and the boundary segment 
   % is oriented correctly.  Otherwise, switch the order of the boundary
   % segment node numbers.  The test points whose corresponding element
   % number, as returned by FINDELE, is NaN is a test point from a boundary
   % segment in the reverse order.
   %
   j=findelemex5(testx,testy,AR,A,B,T);
   
   irev=find(isnan(j));
   temp=bnd(irev,1);
   bnd(irev,1)=bnd(irev,2);
   bnd(irev,2)=temp;

   % Resort Col#1 of bnd and permute Col#2, 
   % just to make it easier to debug
   %
   ib=bnd(:,1);jb=bnd(:,2);
   [ib,iperm]=sort(ib);
   jb=jb(iperm);
   bnd=[ib(:) jb(:)];

return

%
%        Brian O. Blanton
%        Department of Marine Sciences
%        12-7 Venable Hall
%        CB# 3300
%        University of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        brian_blanton@unc.edu
%
%        Summer 2000
%
