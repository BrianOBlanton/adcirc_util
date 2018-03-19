function [LO,LA]=adcirc_wind_rg(temp);
% Call as: [LO,LA]=adcirc_wind_rg([NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC]);

NWLAT =temp(1);
NWLON =temp(2);
WLATMAX =temp(3);
WLONMIN =temp(4);
WLATINC =temp(5);
WLONINC=temp(6);

lo=WLONMIN+(0:NWLON-1)*WLONINC;
la=(WLATMAX-(0:NWLAT-1)*WLATINC)';

LO=repmat(lo,[NWLAT 1]);
LA=repmat(la,[1 NWLON]);

figure
plotcoast('worldcoast')
plotcoast('states')
line(LO',LA','Color','r','LineWidth',.5)
line(LO,LA,'Color','r','LineWidth',.5)  
axeq
axis([-100 -58 5 50])


