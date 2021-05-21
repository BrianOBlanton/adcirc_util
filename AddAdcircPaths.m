function AddAdcircPaths(ADCPTH)

addpath([ADCPTH '/basics'])
addpath([ADCPTH '/owi'])
addpath([ADCPTH '/mex'])
addpath([ADCPTH '/nc'])
addpath([ADCPTH '/viz'])
addpath([ADCPTH '/grid'])
addpath([ADCPTH '/IO'])
addpath([ADCPTH '/misc'])
addpath([ADCPTH '/tides'])   
addpath([ADCPTH '/owi'])   
addpath([ADCPTH '/adcircAnimator'])   

% add javapath for jts
javaaddpath([ADCPTH '/java/jts-1.9.jar']);
