function setup_adcirc_util
global ADCIRC

addpath([ADCIRC '/basics'])
addpath([ADCIRC '/owi'])
addpath([ADCIRC '/mex'])
addpath([ADCIRC '/nc'])
addpath([ADCIRC '/viz'])
addpath([ADCIRC '/grid'])
addpath([ADCIRC '/IO'])
addpath([ADCIRC '/misc'])
addpath([ADCIRC '/tides'])   
addpath([ADCIRC '/owi'])   
addpath([ADCIRC '/adcircAnimator'])   

% add javapath for jts
javaaddpath(genpath([ADCIRC '/java/jts-1.9.jar']));
