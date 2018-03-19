function [x,y]=AdcircCppForward(lambda,phi,lambda0,phi0)

r=6378206.4d0;
x=r*(lambda-lambda0).*cos(phi0);
y=phi*r;

end
