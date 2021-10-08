% MakeRegAdcGrid - class to generate a rectangular ADCIRC grid
% from "OWI" grid parameters {'iLat','iLon', 'dx','dy','SWLat','SWLon'}.
%                 
% Methods: 
%   MakeRegAdcGrid - Instance the class and add coordinates and ADCIRC grid
%                    parts
%   expand_coords  - expand the grid parameters into 2-d lon, lat arrays
%   make_adc_grid  - add element and boundary lists to grid
%   write_grid     - write grid to fort.14 file, currently named "fakegrid.grd"
%   show_grid      - make a simple plot of the simple grid
%
% Usage: 
%   Instance the class and add grid components:  
%       g=MakeRegAdcGrid() returns a default small grid with 
%           iLat: 10, iLon: 10, dx: 10, dy: 10, SWLat: 0, SWLon: 0
%
%       g=MakeRegAdcGrid(struct) returns a grid defined by the grid
%           parameters in the input data struct with fields
%           {'iLat','iLon','dx','dy','SWLat','SWLon'}
%
%       g=MakeRegAdcGrid(owifile) returns a grid defined by the first 
%           grid specification in the owi file owifile.
%           ex: owifile=''/Users/bblanton/Projects/OWI_Winds_5min/fort.221'
%           g=MakeRegAdcGrid(owifile)
%           
%       g=MakeRegAdcGrid([iLat  iLon dx dy  SWLat SWLon]) returns a grid 
%           defined by the 6 grid parameters.  Order is assumed
%           to be correct. 
%
%    Once a grid is "built", write it to disk with 
%    	g.write_grid
% 
%    Viz the simple grid with: 
%		g.show_grid
%

classdef MakeRegAdcGrid

    properties
        file='none'
        iLat = 10
        iLon = 10
        dx = 10
        dy = 10
        SWLat = 0
        SWLon = 0
        Lon
        Lat
        Grid
    end

    methods

        %class init with grid parameters
        function obj=MakeRegAdcGrid(arg)
            
            if ~exist('arg','var')
                
                disp('returning default grid.')
                
            elseif isstruct(arg)
                
                % assume a struct with grid parameters
                needed_fields={'iLat','iLon','dx','dy','SWLat','SWLon'};
                if ~all(isfield(arg,needed_fields))
                    disp('missing data in input struct')
                    return
                end
                obj.iLat=arg.iLat;
                obj.iLon=arg.iLon;
                obj.dx=arg.dx;
                obj.dy=arg.dy;
                obj.SWLat=arg.SWLat;
                obj.SWLon=arg.SWLon;
                
            elseif ischar(arg)
                
                % assume a filename of an owi file, get the grid header
                % line
                if ~exist(arg,'file')
                    error('input file not found. Terminal.')
                end
                fmtstring='iLat=%diLong=%dDX=%fDY=%fSWLat=%fSWLon=%fDT=%4d%2d%2d%2d%2d';
                fid=fopen(arg,'r');
                headerline=fgets(fid);
                l=fgets(fid);
                fclose(fid);
                A=sscanf(l,fmtstring);
                obj.iLat=A(1);
                obj.iLon=A(2);
                obj.dx=A(3);
                obj.dy=A(4);
                obj.SWLat=A(5);
                obj.SWLon=A(6);
                obj.file=arg;
                
            elseif length(arg)==6
                
                % assume order
                obj.iLat=arg(1);
                obj.iLon=arg(2);
                obj.dx=arg(3);
                obj.dy=arg(4);
                obj.SWLat=arg(5);
                obj.SWLon=arg(6);
                
            else
                
                error('unrecognized input to MakeRegAdcGrid.')
                
            end
            
            obj=obj.expand_coords();
            obj=obj.make_adc_grid();
            
        end

        function obj=expand_coords(obj)
            
            XGrid=obj.SWLon+(0:obj.iLon-1)*obj.dx;
            YGrid=obj.SWLat+(0:obj.iLat-1)*obj.dy;
            [obj.Lon,obj.Lat]=meshgrid(XGrid,YGrid);
            
        end
        
        function obj=make_adc_grid(obj)
            
            [nx,ny]=size(obj.Lon);
            g.name='fakegrid';
            g.x=obj.Lon(:);
            g.y=obj.Lat(:);
            g.z=nan(size(g.y));
            g.e=elgen(ny,nx);
            g.bnd=detbndy(g.e);
            obj.Grid=g;
            
        end
        
        function write_grid(obj, name)
            
            if exist('name','var')
                obj.Grid.name=name;
            end
            
            fgs2fort14(obj.Grid)
            
        end
        
        function show_grid(obj)
            
            figure
            plotbnd(obj.Grid)
            drawelems(obj.Grid)
            line(obj.Grid.x,obj.Grid.y,'marker','.','linestyle','none')
            line(obj.Grid.x',obj.Grid.y','marker','.','linestyle','none')
            axis('equal')
            title(obj.Grid.name)

        end
        
    end
    
end

