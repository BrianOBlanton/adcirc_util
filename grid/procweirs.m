function procweirs(fgs)

%
% Generate a weir data table for RA model
%

%
% IPET RA work
% Brian Blanton, SAIC
% 1 Sep 2006
%

fid=fopen('weir.dat','w');
fprintf(fid,'%s\n',fgs.name);
fprintf(fid,'   BND#   SEG#    LON1       LAT1       LON2        LAT2       LEN      H1      H2       Z1     Z2\n');
fprintf(fid,'                  deg         deg       deg         deg         ft      ft      ft       ft     ft\n');

LO0=-90.1884;
LA0=29.9153;

for i=1:fgs.nland
   switch fgs.ibtype(i)
   case 24
	   for j=1:length(fgs.ln{i})
		  
		  n1=fgs.ln{i}(j,1);
		  n2=fgs.ln{i}(j,2);
		  lo1=fgs.x([n1]);
		  lo2=fgs.x([n2]);
		  la1=fgs.y([n1]);
		  la2=fgs.y([n2]);
		
		  [x1,y1]=convll2m(lo1,la1,LO0,LA0);
		  [x2,y2]=convll2m(lo2,la2,LO0,LA0);

		  h1=fgs.weirheights{i}(j);
		  h2=fgs.weirheights{i}(j);
		  z1=fgs.z(n1);  
		  z2=fgs.z(n2);  
		  % 1 meter = 3.280839895 feet
		  len=sqrt((x2-x1).^2 + (y2-y1).^2)*3.280839895;
		  len2=sqrt((lo2-lo1).^2 + (lo2-lo1).^2);
		  out=[i j lo1 la1 lo2 la2 len [h1 h2 z1 z2]*3.280839895];

		  fprintf(fid,'%6d %6d %10.4f %10.4f %10.4f %10.4f %10.2f %7.2f %7.2f %7.2f %7.2f\n',out);

	   end
	end
end

fclose(fid);


