function [x,y]=AdcircCppForward(lambda,phi,lambda0,phi0)
% AdcircCppForward ADCIRC's Platte Carre projection
%
% Inputs:
%   lambda  : longitude in degrees
%   phi     : latitude in degrees
%   lambda0 : longitude at center of projection (degrees)
%   phi0    : latitude at center of projection (degrees)
%
% Outputs:
%   x       : projected longitude (m)
%   y       : projected latitude (m)
%
%
% Call as: [x,y]=AdcircCppForward(lambda,phi,lambda0,phi0);

% if range(lambda)<1
%     fprintf('Input coordinates might be in radians.\nIf you get bizarre results, this might be why.\n')
% end

r=6378206.4d0;
x=r*(lambda-lambda0)*pi/180.*cos(phi0*pi/180);
y=phi*pi/180*r;

