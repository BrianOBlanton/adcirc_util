function showmax(g,D)
% showmax(g,data)

figure
%drawelems(g,'Color',[1 1 1]*.7,'Linewidth',.25);
plotbnd(g,'LineWidth',2);
lcontour(g,'z',0,'Color','k','LineWidth',.2);
axis('equal')

colormesh2d(g,D)
colorbar

[~,b]=max(D);

line(g.x(b),g.y(b),10,'Marker','*','Color','g','MarkerSize',20)