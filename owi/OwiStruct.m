function OwiStruct
%  A structure with the following fields: 
%
%           Basin: [1x1 struct]
%                   Pre:  {1xNT cell}
%                   WinU: {1xNT cell}
%                   WinV: {1xNT cell}
%                   Time: [1xNT double]
%           Region: [1x1 struct]
%                   Pre:  {1xNT cell}
%                   WinU: {1xNT cell}
%                   WinV: {1xNT cell}
%                   Time: [1xNT double]
%
%          where NT is the number of times in the OWI struct. Each cell 
%          value is the size of X (and Y). The Region field will exist 
%          only if the OWI struct contains a Region field.  The time
%          lengths of the Basin and Region fields must be the same.  
%
% To convert the cell array into a double array, do: 
%
% BasinPre=reshape([D.Basin.Pre{:}],[size(X) length(OWI.Basin.time)]);
%
% If X is 2-D, then it may ne necessary to "squeeze" a dimension out of 
% BasinPre.  E.g., to plot the timeseries for the first point, do:
%
% plot(OWI.Basin.time,squeeze(BasinPre(1,1,:)))
%
str={
'A structure with the following fields: ',
' ',
'          Basin: [1x1 struct]',
'                  Pre:  {1xNT cell}',
'                  WinU: {1xNT cell}',
'                  WinV: {1xNT cell}',
'                  Time: [1xNT double]',
'          Region: [1x1 struct]',
'                  Pre:  {1xNT cell}',
'                  WinU: {1xNT cell}',
'                  WinV: {1xNT cell}',
'                  Time: [1xNT double]',
' ',
'         where NT is the number of times in the OWI struct. Each cell ',
'         value is the size of X (and Y). The Region field will exist ',
'         only if the OWI struct contains a Region field.  The time',
'         lengths of the Basin and Region fields must be the same.  ',
' ',
'To convert the cell array into a double array, do: ',
' ',
'BasinPre=reshape([D.Basin.Pre{:}],[size(X) length(OWI.Basin.time)]);',
' ',
'If X is 2-D, then it may ne necessary to "squeeze" a dimension out of',
'BasinPre.  E.g., to plot the timeseries for the first point, do:',
' ',
'plot(OWI.Basin.time,squeeze(BasinPre(1,1,:)))'
};
sprintf('%s\n',str{:})