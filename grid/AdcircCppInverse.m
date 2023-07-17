function [lambda,phi] = AdcircCppInverse(x,y,lambda0, phi0)
% AdcircCppInverse Inverse of ADCIRC's Platte Carre projection
%
% Inputs:
%   x       : projected longitude (m)
%   y       : projected latitude (m)
%   lambda0 : longitude at center of projection (degrees)
%   phi0    : latitude at center of projection (degrees)
%
% Outputs:
%   lambda  : longitude in degrees
%   phi     : latitude in degrees
%
% Call as: [lambda,phi] = AdcircCppInverse(x,y,lambda0,phi0)

r=6378206.4d0;
alpha=cos(phi0*pi/180);
lambda=lambda0+180/pi*(x./(r*alpha));
phi=180/pi*y/r;

