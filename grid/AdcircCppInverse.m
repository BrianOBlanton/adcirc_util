function [lambda,phi] = AdcircCppInverse( x,y,lambda0, phi0)
% AdcircCppInverse Summary of this function goes here
%   Detailed explanation goes here


r=6378206.4d0;
lambda=rlambda0+x./(r.*cos(phi0));
phi=y/r;



end

