classdef nws13

    properties
        url='None';
        % grids contains:
        % lon, lat, U10, V10, PSFC, time
        grids
        groups
        ax
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

                i=find(obj.grids(j).time==t);

                lon=squeeze(obj.grids(j).lon(i,:,:));
                lat=squeeze(obj.grids(j).lat(i,:,:));
                psfc=squeeze(obj.grids(j).PSFC(i,:,:));

                pcolor(lon,lat,psfc)

                xb=[lon(1,:)'
                    lon(:,end)
                    lon(end,:)'
                    flipud(lon(:,1))];
                yb=[lat(1,:)'
                    lat(:,end)
                    lat(end,:)'
                    flipud(lat(:,1))];

                line(xb,yb,Color='k')

            end

        end
            
        function obj=drawMain(obj,t)

            i=find(obj.grids(1).time==t);

            hold on
            pcolor(obj.grids(1).lon(:,:), obj.grids(1).lat(:,:), ...
                   squeeze(obj.grids(1).PSFC(i,:,:)))
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
