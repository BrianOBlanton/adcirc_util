function [lambda,phi] = AdcircCppInverse( x,y,lambda0, phi0)
% AdcircCppInverse Inverse of ADCIRC's Platte Carre projection
% 
% Call as: [lambda,phi] = AdcircCppInverse(x,y,lambda0, phi0)


r=6378206.4d0;
lambda=lambda0+x./(r.*cos(phi0));
phi=y/r;



end

