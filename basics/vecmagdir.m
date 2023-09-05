function [vm,vd]=vecmagdir(u,v)
%VECMAGDIR compute magnitude and direction of vector
% VECMAGDIR computes the magnitude and direction of a vector
% vector field v.
%
%  Inputs: u,v - vector field
% Outputs: vm - magnitude
%          vd - direction in DEGREES!!!
%
% Call as: [vm,vd]=vecmagdir(u,v);

% BOB, Summer 2000


if nargin==0
    disp('[vm,vd]=vecmagdir(u,v);');
    return
end

vm=sqrt(u.*u+v.*v);
vd=atan2(v,u)*180/pi;
