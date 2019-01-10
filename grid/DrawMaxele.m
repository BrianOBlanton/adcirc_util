function DrawMaxele(nc)
% DrawMaxele(<maxele.63.nc>)

if nargin==0
    disp('Dout=DrawMaxele(<maxele.63.nc>)')
    return
end

figure

if ~isa(nc,'ncgeodataset')
   %assume its a filename to load
   nc=ncgeodataset(nc);
end

FGS=ExtractGrid(nc);
r=range(FGS.x);
zz=nc{'zeta_max'}(:);

[~,b]=max(zz);

ax=[FGS.x(b)-r/20 FGS.x(b)+r/20 FGS.y(b)-r/20 FGS.y(b)+r/20];

colormesh2d(FGS,zz);
plotbnd(FGS,'Color','k','LineWidth',2)
lcontour(FGS,'z',0,'LineWidth',.25,'Color','k');
axis('equal')

caxis([0 ceil(3*nanstd(zz))])
colormap(jet(20))
hcb=colorbar;
hcb.Label.String='m [MSL]';
grid on
grid minor