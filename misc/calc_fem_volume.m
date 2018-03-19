function v=calc_fem_volume(fem_grid_struct,zeta)
%CALC_FEM_VOLUME - ESTIMATE!! the volume in a FEM domain
% CALC_FEM_VOLUME estimates the total volume on a FEM domain
%
%  Inputs: fem_grid_struct
%          zeta - surface elevation
%  Output: v - total volume estimate
%
% Call as: v=calc_fem_volume(fem_grid_struct,zeta);
%


% Area times average depth over element

d=-fem_grid_struct.z-zeta;

d3=mean(d(fem_grid_struct.e)')';

v=sum(abs(d3).*fem_grid_struct.ar);
