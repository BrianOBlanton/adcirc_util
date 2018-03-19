function Dout=DrawMaxele(FGS,Din)
% Dout=DrawMaxele(FGS,Din)

if nargin==0
    disp('Dout=DrawMaxele(FGS,Din)')
    return
end

%figure

if ~exist('Din')
   Din=read_adcirc_fort('FileName','maxele.63','Compact',0)
   Dout=Din;
end

colormesh2d(FGS,Din.zeta);
plotbnd(FGS)
lcontour(FGS,'z',0,'LineWidth',.25,'Color','k');
axis('equal')

colorbar

grid