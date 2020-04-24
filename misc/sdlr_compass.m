function h=sdlr_compass(f13,n,varargin)
% function h=sdlr_compass(f13,n,varargin)
% Inputs: f13 - complete fort.13 file struct
%         n   - node to plot
% Outputs: handles to vectors drawn.


theta=(15:30:370)*pi/180;
if n>f13.nn
    error('Input node number (%d) exceeds size of grid (%d).',n,f13.nn)
end

isdrl=find(strcmp('surface_directional_effective_roughness_length',f13.att_names));

uv=f13.atts{isdrl}(n,:).*exp(1i*theta);

h=compass(uv,varargin{:});

