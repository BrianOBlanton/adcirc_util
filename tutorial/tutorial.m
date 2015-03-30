%% Handling ADCIRC files in MATLAB
% AdcMat - a MATLAB toolbox for ADCIRC

% Brian Blanton
% Renaissance Computing Institute
% The University of North Carolina at Chapel Hill
% brian_blanton@renci.org

%% Loading an ADCIRC grid
% fgs=grd_to_opnml(adcirc_grid_name); 
%
% where adcirc_grid_name is a .grd file or fort.14, and fgs is the output 
% variable containing the ADCIRC grid. The data returned by GRD_TO_OPNML is
% a MATLAB structure called a fem_grid_struct. Type "help fem_grid_struct" 
% for more information. Most ADCIRC-related plotting functions require a 
% fem_grid_struct as input. 

fgs=grd_to_opnml('ec95d.grd');

%% Drawing the ADCIRC grid
% he=drawelems(fem_grid_struct,pn1,pv1,...); 
% hb=plotbnd(fem_grid_struct,pn1,pv1,...); 
%	

figure
he=drawelems(fgs,'Color','b','LineWidth',.1);
hb=plotbnd(fgs,'Color','k','LineWidth',2);
axis('equal')
snapnow

% The figure can then be zoomed-in, scaled, annotated, etc, as usual.  For example: 
axis([-78.1  -72.8   33.3   37.1])


%% HTML Markup Example
% This is a table:
%
% <html>
% <table border=1><tr><td>one</td><td>two</td></tr>
% <tr><td>three</td><td>four</td></tr></table>
% </html>
%
