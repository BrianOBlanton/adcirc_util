function OwiStruct
%  A structure with the following 2 fields, each of which as more fields:
%
%           Basin: [1x1 struct]
%                   time: [1xNT double]   % Time vector
%                   iLat: [1xNT double]   % Number of latitude points
%                  iLong: [1xNT double]   % Number of latitude points
%                     DX: [1xNT double]   % Longitude delta
%                     DY: [1xNT double]   % Latitude delta
%                  SWLat: [1xNT double]   % Latitude of southwest corner of grid
%                  SWLon: [1xNT double]   % Longitude of southwest corner of grid
%                  XGrid: {1xNT cell}     % Longitude of grid (for each time level)
%                  YGrid: {1xNT cell}     % Latitude of grid (for each time level)
%                    Pre: {1xNT cell}     % Pressure on grid (for each time level)
%                   WinU: {1xNT cell}     % East/west wind speed on grid (for each time level)
%                   WinV: {1xNT cell}     % North/south wind speed on grid (for each time level)
%                 
%           Region: [1x1 struct]
%                   same as above
%                 
%          where NT is the number of times in the OWI fields. Each cell 
%          value is the size of X (and Y). The Region field will exist 
%          only if the OWI struct contains a Region field.  The time
%          lengths of the Basin and Region fields must be the same.  
%
% To convert the cell array into a double array, do: 
%
% BasinPre=reshape([D.Basin.Pre{:}],[size(X) length(OWI.Basin.time)]);
%
% If X is 2-D, then it may be necessary to "squeeze" a dimension out of 
% BasinPre.  E.g., to plot the timeseries for the first point, do:
%
% plot(OWI.Basin.time,squeeze(BasinPre(1,1,:)))
%
