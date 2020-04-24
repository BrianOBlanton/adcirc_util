function h=sdelr_compass(f13,n,varargin)
% function h=sdelr_compass(f13,n,varargin)
%
% sdelr_compass draws the roughness length directions from a fort.13 file
% (as read in by read_adcirc_fort13) at a specific list of nodes.
%
% This nodal attribute is 
% directional, and the twelve values represent the roughness 
% lengths ?seen? by winds blowing from twelve different compass 
% directions at each node. The orientation of the twelve values 
% follows the trigonometric convention, that is, zero degrees 
% represents due east, and the values proceed counter clockwise. 
% In other words, the first value at a node is applied to winds 
% blowing from west to east, the second value applies to winds 
% blowing East-Northeast, etc.
%
% Inputs: f13 - complete fort.13 file struct
%         n   - nodes to plot
% Outputs: handles to vectors drawn.


theta=(15:30:370)*pi/180;
if n>f13.nn
    error('Input node number (%d) exceeds size of grid (%d).',n,f13.nn)
end

isdrl=find(strcmp('surface_directional_effective_roughness_length',f13.att_names));

uv=f13.atts{isdrl}(n,:).*exp(1i*theta);

h=compass(uv,varargin{:});

