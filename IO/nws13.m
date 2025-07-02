classdef nws13

    properties
        url='None';
        % grids contains:
        % lon, lat, U10, V10, PSFC, time
        grids
        groups
        ax
    end

    properties (Constant)
        nest_vec_stride=5;
        main_vec_stride=1;
    end

    methods

        %class init
        function obj=nws13(url)
            if ~exist(url)
                error('Input url DNE.')
            end
            %obj=struct();
            obj.url=url;
            obj=load(obj);
        end

        function obj=load(obj)

            nc=ncgeodataset(obj.url);

            groupnames=split(nc.attribute{'group_order'});
            ngroups=length(groupnames);
            varnames={'lon' 'lat' 'PSFC' 'U10' 'V10' 'time'};
        
            obj.groups=groupnames;

            for i=1:ngroups
                obj.grids(i).group=groupnames{i};
                for j=1:length(varnames)
                    vstr=sprintf('%s/%s',groupnames{i},varnames{j});
                    obj.grids(i).(varnames{j})=nc{vstr};
                    if strcmp(varnames{j},'time')

                        ncc=nc{vstr};
                        units=split(ncc.attribute('units'));
                        stime=datetime(units{3});
                        times=ncc.data(:);

                        switch lower(units{1})
                            case 'minutes'
                                obj.grids(i).time=stime+minutes(times);
                            case 'seconds'
                                obj.grids(i).time=stime+seconds(times);
                            otherwise
                                error('uncoded time units: %s',units{1})
                        end
                        
                    end
                end

                obj.ax=[obj.grids(1).lon(1,1) obj.grids(1).lon(1,end) obj.grids(1).lat(end,end) obj.grids(1).lat(1,1)];
            
            end
        end

        function obj=drawSnap(obj,t)

            obj=drawMain(obj,t);
            obj=drawNests(obj,t);
            addStuff(obj)
            title(string(t),'FontSize',30)

        end

        function obj=drawNests(obj,t)

            for j=2:length(obj.grids)

                i=find(obj.grids(j).time==t)

                if ~isempty(i)

                    lon= squeeze(obj.grids(j).lon(i,1:obj.nest_vec_stride:end,1:obj.nest_vec_stride:end));
                    lat= squeeze(obj.grids(j).lat(i,1:obj.nest_vec_stride:end,1:obj.nest_vec_stride:end)); 
                    psl=squeeze(obj.grids(j).PSFC(i,1:obj.nest_vec_stride:end,1:obj.nest_vec_stride:end)); 
                    u10= squeeze(obj.grids(j).U10(i,1:obj.nest_vec_stride:end,1:obj.nest_vec_stride:end)); 
                    v10= squeeze(obj.grids(j).V10(i,1:obj.nest_vec_stride:end,1:obj.nest_vec_stride:end)); 
                  
                    pcolor(lon,lat,psl)
                    hold on
                    hv=vecplot(lon,lat,u10,v10,'ScaleFac',10,'ScaleLabel','no scale');
%                    hv=quiver(lon,lat,u10,v10,10);

                    hv.Color='k';

                    xb=[lon(1,:)'
                        lon(:,end)
                        lon(end,:)'
                        flipud(lon(:,1))];
                    yb=[lat(1,:)'
                        lat(:,end)
                        lat(end,:)'
                        flipud(lat(:,1))];
    
                    line(xb,yb,Color='k')
                    text(xb(1),yb(1),string(j))

                else

                    fprintf('Nothing to draw for rank %d at time %s.\n',j,string(t))

                end
            end

        end
            
        function obj=drawMain(obj,t)

            i=find(obj.grids(1).time==t);

            hold on
            pcolor(obj.grids(1).lon(:,:), ...
                   obj.grids(1).lat(:,:), ...
                   squeeze(obj.grids(1).PSFC(i,:,:)))

            % hv=quiver(obj.grids(1).lon(1:obj.main_vec_stride:end,1:obj.main_vec_stride:end),...
            %           obj.grids(1).lat(1:obj.main_vec_stride:end,1:obj.main_vec_stride:end),...
            %           squeeze(obj.grids(1).U10(i,1:obj.main_vec_stride:end,1:obj.main_vec_stride:end)),...
            %           squeeze(obj.grids(1).V10(i,1:obj.main_vec_stride:end,1:obj.main_vec_stride:end)),10);
            hv=vecplot(obj.grids(1).lon(1:obj.main_vec_stride:end,1:obj.main_vec_stride:end),...
                      obj.grids(1).lat(1:obj.main_vec_stride:end,1:obj.main_vec_stride:end),...
                      squeeze(obj.grids(1).U10(i,1:obj.main_vec_stride:end,1:obj.main_vec_stride:end)),...
                      squeeze(obj.grids(1).V10(i,1:obj.main_vec_stride:end,1:obj.main_vec_stride:end)),...
                      'ScaleFac',100,'ScaleLabel','no scale');


            hv.Color='k';


        end

        function obj=addStuff(obj)
            colormap(jet(20))
            colorbar
            shading interp
            ax=gca;
            ax.Layer='top';
            ax.GridLineStyle='-';
            ax.MinorGridLineStyle=':';
            ax.GridAlpha=.5;
            ax.MinorGridAlpha=.5;
            grid on
            grid minor
            ax.Clipping='on';
            axis('equal')
            axx=axis;
            drawgshhg;
            axis(axx)
        end
    end
end
