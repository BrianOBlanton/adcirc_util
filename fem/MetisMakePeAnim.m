
pes=[64 65 66 67 69 70 71 76 77 79];

npes=length(pes);

ShowFort69=false;
ShowFort63=true;
ShowFort64=false;
ShowFort74=false;
LabelPEs=false;
DepthContours=[-2 -1 1];
Stride=1;
Start=1;
End=[];

%% load master grid
FGS=grd_to_opnml;

%%
for i=pes
   com=sprintf('g%04d=grd_to_opnml(''PE%04d/fort.14'');',i,i);
   eval(com)
   com=sprintf('g%04d.mx=mean(g%04d.x);',i,i);   
   eval(com)
   com=sprintf('g%04d.my=mean(g%04d.y);',i,i);
   eval(com)
   if ShowFort63
      disp(sprintf('   Loading PE%04d/fort.63',i))
      com=sprintf('D63_%04d=read_adcirc_fort(''FileName'',''PE%04d/fort.63'',''FortUnit'',''63'',''Compact'',0);',i,i);
      eval(com)
   end
   if ShowFort64
      disp(sprintf('   Loading PE%04d/fort.64',i))
      com=sprintf('D64_%04d=read_adcirc_fort(''FileName'',''PE%04d/fort.64'',''FortUnit'',''64'',''Compact'',0);',i,i);
      eval(com)
   end
   if ShowFort69
      disp(sprintf('   Loading PE%04d/fort.69',i))
      com=sprintf('D69_%04d=read_adcirc_fort(''FileName'',''PE%04d/fort.69'',''FortUnit'',''63'',''Compact'',0);',i,i);
      eval(com)   
   end
   if ShowFort74
      disp(sprintf('   Loading PE%04d/fort.74',i))
      com=sprintf('D74_%04d=read_adcirc_fort(''FileName'',''PE%04d/fort.74'',''FortUnit'',''74'',''Compact'',0);',i,i);
      eval(com) 
   end
end


%%

ax=[  -78.7162  -77.1622   33.7201   34.5085];
ax=[  -78.4738  -77.5057   33.8553   34.3465];
ax=[  -78.1900  -77.7397   34.0945   34.3488];
ax=[  -77.9924  -77.9471   34.2382   34.2612];


if ShowFort63
   FirstPEFort=eval(sprintf('D63_%04d',pes(1)));
elseif ShowFort69
   FirstPEFort=eval(sprintf('D69_%04d',pes(1)));
end

FirstPEGrid=eval(sprintf('g%04d',pes(1)));

if ~exist('End') | isempty(End)
   End=length(FirstPEFort.t);
end
End=min(End,length(find(isfinite(FirstPEFort.t))));

Start=2560;
Stride=1;
End=3000;

%% set up figure

figure

%colormap(tmap)
axeq
axis(ax)
cm=max(d1);
%set(gca,'CLim',cm*[-1 1])
set(gca,'CLim',[-2 1])
colorbar
colormap(jet(24))
hel=drawelems(FGS,'Color','k','Linewidth',.25);
set_height(hel,1);
hz0=lcontour(FGS,'z',0,'Color','k','LineWidth',1);
set_height(hz0,1);

%%
clear h63 h64 h69 h74
d1=NaN*ones(length(pes),1);

for idx=Start:Stride:End
    
   %if isnan(FirstPEFort.t(idx))
   %   break
   %end
   
   % delete prev graphics objects
   if exist('h63'),delete(h63);h63=NaN*ones(npes,1);,end
   if exist('h64'),delete(h64);h64=NaN*ones(npes,1);,end
   if exist('h69'),delete(h69);h69=NaN*ones(npes,1);,end   
   if exist('h74'),delete(h74);h74=NaN*ones(npes,1);,end   
   
   k=1;
   
   for j=pes
       
      if ShowFort63
         com=sprintf('h63(%d)=colormesh2d(g%04d,D63_%04d.zeta(:,idx));plotbnd(g%04d);',k,j,j,j);
         eval(com)   
         com=sprintf('d1(%d)=max(abs(D63_%04d.zeta(:,idx)));',k,j);
         eval(com)
      end
      
      if ShowFort64
        com=sprintf('h64(%d)=vecplot(g%04d.x,g%04d.y,D64_%04d.ubar(:,idx),D64_%04d.vbar(:,idx),''ScaleLabel'',''no scale'',''ScaleFac'',1,''Stride'',1,''Color'',''k'');',k,j,j,j,j);
         eval(com)
      end

      if ShowFort69
         com=sprintf('h69(%d)=colormesh2d(g%04d,D69_%04d.zeta(:,idx));plotbnd(g%04d);',k,j,j,j);
         eval(com)
         
         com=sprintf('d1(%d)=max(abs(D69_%04d.zeta(:,idx)));',k,j);
         eval(com)
      end

      if ShowFort74
         com=sprintf('h74(%d)=vecplot(g%04d.x,g%04d.y,D74_%04d.uwind(:,idx),D74_%04d.vwind(:,idx),''ScaleLabel'',''no scale'',''ScaleFac'',50,''Stride'',50,''Color'',''k'');',k,j,j,j,j);
         eval(com)
      end
      
     
      
      %com=sprintf('lcontour(g%04d,''z'',DepthContours,''LineWidth'',.25);',j);
      %eval(com)
      %com=sprintf('lcontour(g%04d,''z'',0,''Color'',''k'',''LineWidth'',.25);',j);
      %eval(com)
      if LabelPEs
         com=sprintf('text(g%04d.mx,g%04d.my,''%04d'',''FontSize'',8,''FontWeight'',''normal'');',j,j,j);
         eval(com)
      end
      
      k=k+1;
   
   end
         
   d2=max(d1);
   %tstr{1}=sprintf('t=%s',datestr(FirstPEFort.t(idx)/86400,0));
   tstr{1}=sprintf('t=%d days  %s   Max=%6.2f',floor(FirstPEFort.t(idx)/86400), datestr(FirstPEFort.t(idx)/86400,13),d2);
   tstr{2}=sprintf('frame=%04d  :  iter=%d',idx,FirstPEFort.iter(idx));

   title(tstr)

   %stamp_right(int2str(idx));
   drawnow
   fname=sprintf('frame.%04d',idx);
   print('-dpng','-r200',[fname '.png']) 
   %ConvertImage([fname '.png'],'gif');

end
return

%%
figure
for idx=1:24
   clf
   for j=pes
      com=sprintf('colormesh2d(g0%d,D69_0%d.zeta(:,idx));plotbnd(g0%d);',j,j,j);
      eval(com)
   end
   %lcontour(g2,'z',[0],'Color','r','Linewidth',.25)
   title(sprintf('Fort.69 @ t=%.1f',D69_0424.t(idx)))
   axeq
   axis([  -77.7061  -77.5295   34.3226   34.4464])
   colorbar
   drawnow
   fname=sprintf('frame.69.%02d',idx);
   print('-dpng','-r200',[fname '.png']) 
   eval(sprintf('!/usr/local/bin/convert %s.png %s.gif',fname,fname));

end

