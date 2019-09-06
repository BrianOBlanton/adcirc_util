function fgs2fort14(fgs)
%FGS2FORT14 output a fort.14 file from a fem_grid_struct
% Call as: fgs2fort14(fgs)

nn=size(fgs.x,1);
ne=size(fgs.e,1);

fid=fopen([fgs.name '.grd'],'w');

fprintf(fid,'%s\n',fgs.name);
fprintf(fid,'%d %d\n',size(fgs.e,1),size(fgs.x,1));

% node list 
out=[(1:nn)' fgs.x fgs.y fgs.z];
fprintf(fid,'%10d %f %f %10.4f\n',out');

% element list 
out=[(1:ne)' 3*ones(size(fgs.e,1),1) fgs.e];
fprintf(fid,'%10d %2d %10d %10d %10d\n',out');

% total number of open boundary nodes
if isfield(fgs,'nopenboundaries')
   fprintf(fid,'%-20d ! Number of open boundaries\n',fgs.nopenboundaries);  % nope   
   for i=1:fgs.nopenboundaries
       fprintf(fid,'%-20d ! Total number of open boundary nodes\n',fgs.nopenboundarynodes{i});  % neta, total number of elevation nodes
       %fprintf(fid,'%-20d Total number of open boundary nodes\n',fgs.elevation);  % neta, total number of elevation nodes
       fprintf(fid,'%-20d\n',fgs.nopenboundarynodes{i});
       fprintf(fid,'%-20d \n',fgs.ob{i});
   end
else
   fprintf(fid,'%-20d ! Number of open boundaries\n',0);  % nope   
   fprintf(fid,'%-20d ! Total number of open boundary nodes\n',0);  % neta, total number of elevation nodes
end

% total number of no-normal-flow
if isfield(fgs,'nland')
      fprintf(fid,'%-20d ! Number of land boundaries\n',fgs.nland);  % default nbou   
      fprintf(fid,'%-20d ! Total number of flux boundary nodes\n',fgs.nlandnodestotal);  % default nbou 
      for i=1:fgs.nland
           fprintf(fid,'%-20d %2d  !Nodes for land bnd %d\n',fgs.nlandnodes(i),fgs.ibtype(i),i); 
           fprintf(fid,'%-10d\n',fgs.ln{i});
      end
else
    fprintf(fid,'%-20d ! Number of land boundaries\n',1);  % default nbou
    fprintf(fid,'%-20d ! Total number of flux boundary nodes\n',size(fgs.bnd,1));  % default nbou
    fprintf(fid,'%-10d%-10d ! Number of nodes for boundary 1\n',size(fgs.bnd,1),0);  % first and only lnd bndy
    out=unique(fgs.bnd);
    fprintf(fid,'%-10d\n',out);
end

fclose(fid);

return
