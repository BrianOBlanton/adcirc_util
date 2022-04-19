function AddAdcircPaths(ADCPTH)
addpath(genpath([ADCPTH]))
return

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
addpath([ADCPTH '/extern'])   

% add javapath for jts
javaaddpath([ADCPTH '/java/jts-1.9.jar']);
