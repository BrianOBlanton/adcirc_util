function AddAdcircPaths(ADCPTH)
%addpath(genpath([ADCPTH]))
%return

addpath([ADCPTH '/basics'],'-end')
addpath([ADCPTH '/owi'],'-end')
addpath([ADCPTH '/mex'],'-end')
addpath([ADCPTH '/nc'],'-end')
addpath([ADCPTH '/viz'],'-end')
addpath([ADCPTH '/grid'],'-end')
addpath([ADCPTH '/IO'],'-end')
addpath([ADCPTH '/misc'],'-end')
addpath([ADCPTH '/tides'],'-end') 
addpath([ADCPTH '/owi'],'-end')   
addpath([ADCPTH '/adcircAnimator'],'-end')   
addpath([ADCPTH '/extern'],'-end')   

% add javapath for jts
javaaddpath([ADCPTH '/java/jts-1.9.jar']);
