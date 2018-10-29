function D=drog2ddtnc(t1d,t2d,dt,idt,xi,yi,Z,V,options)
%function D=drog2ddtnc(TheGrid,t1d,t2d,dt,idt,xi,yi,V,options)
%DROG2DDT track drogues in a 2-D FEM domain, time-stepping version
% DROG2DDT tracks particles through a discrete squence of 2-D velocity
% fields, vertically averaged for example.  The integrator is a 2nd order
% Runge-Kutta (mid-point) method.  Pass empty ([]) for default parameter
% values.
%
% Inputs: t1,t2    - integration end-points; both t1 & t2 must lie within
%                    the min and max time in the velocity sequence
%                    (obviously). The time must be in MATLAB's datenum or
%                    datetime format. Defaults= beginning and end of input
%                    velocity field.
%         dt       - integration time step; this need NOT match the
%                    velocity sequence step. This needs to be a MATLAB
%                    datetime duration, e.g., duration(0,15,0) for a 15
%                    minute integration step size. Default=time interval
%                    in velocity field.
%         idt      - output interval; this is the "frequency" with which
%                    to store computed step.  idt=1 means store each step;
%                    idt=2 means store every other step, etc. Default=1.
%         xi,yi    - initial drogue positions; default=10 randomly chosen
%                    grid node locations
%         V        - ncgeodataset object, e.g.,
%                    V=ncgeodataset('fort.64.nc');  Default: local
%                    fort.64.nc.  Can be a 3d output file (fort.45.nc), in
%                    which case the options field iv must be provided if
%                    not "surface"
%         options  - options structure; valid fields are:
%                       'draw'=true - plot updated locations on figure
%                       'verbose'=true - for runtime diags
%                       'lag'=number - plotting tail length in steps
%                       'integrator'='rk4' - ('rk4'|'rk2') - numerical integrator to use
%                       'conv'=do NOT convert grid to CPP.
%                       'iv'=vertical level to integrate (default=surface)
%
% Outputs: A struct containing:
%         xx,yy    - arrays of drogue trajectories.  The size of the
%                    arrays is the size of xi,yi, plus an extra
%                    dimension for the time history.  This dimension
%                    is the number of iterations used, possibly
%                    subsampled at idt intervals.  Example: if
%                    xi is 10x10, and t1=0,t2=10,dt=1,idt=2, then the
%                    size of xx is 10x10x6.
%
%         tt       - times at which outputs are stored
%         uu,vv    - along-track velocity in the same format as xx,yy
%
% Call as: D=drog2ddtnc(t1,t2,dt,idt,xi,yi,V,options);
%

%        Brian O. Blanton
%        RENCI
%        University of North Carolina
%        Chapel Hill, NC
%
%        brian_blanton@renci.org
%
%        March 2018
%           added rk4 integrator
%           mod'd for ncgeodataset input of velocity
%

if nargin==0
    disp('Call as:  D=drog2ddtnc(t1,t2,dt,idt,xi,yi,V,options)')
    return
end

% if nargin==1 && strcmp(fem_grid_struct,'velhelp')
%    velhelp('Help');
%    return
% end

if nargin<8
    error('Too few arguments.')
end

% % First argument to drog2ddtnc must be a grid structure
% if ~is_valid_struct(fem_grid_struct)
%    error('fem_grid_struct to DROG2DDT NOT valid.')
% end

if ~isa(V,'ncgeodataset')
    error('Velocity field input must be an ncgeodataset.')
end

% Check options structure
if isempty(options)
    options.draw=0;
    options.lag=48;
    options.iv=0;
elseif ~isstruct(options)
    error('Options argument to DROG2DDTNC must be a structure')
else
    % check fields
end

global dim3d
dim3d=0;

% Verify velocity input
if any(ismember({'u-vel3D','v-vel3D'},V.variables))
    fprintf('Input ncgeodataset appears to be a 3-d ADCIRC file.\n');
    if ~all(ismember({'u-vel3D','v-vel3D','time'},V.variables))
        error('Needed variables {''u-vel3D'',''v-vel3D'',''time''} are missing in the ncgeodataset object')
    end
    s=size(V{'u-vel3D'});
    if options.iv==0
        dim3d=s(2);
    else
        dim3d=options.iv;
    end
else % 2-d file
    fprintf('Input ncgeodataset appears to be a 2-d ADCIRC file.\n');
    if ~all(ismember({'u-vel','v-vel','time'},V.variables))
        error('Needed variables {''u-vel'',''v-vel'',''time''} are missing in the ncgeodataset object')
    end
end

time=V.time{'time'};  % extract time to datenum
% convert to datetime
time=datetime(datevec(time));

% extract fem_grid_struct from ncgeodataset
TheGrid=ExtractGrid(V);
ConvToXY=false;
% convert to cart if needed.
if max(TheGrid.x)-min(TheGrid.x)<50
    TheGridll=TheGrid;
    fprintf('   Grid appears to be in Lon/Lat. Converting to CPP.\nIf this is wrong, pass in opts.conv=false.\n')
    ConvToXY=true;
    lo0=mean(TheGrid.x);
    la0=mean(TheGrid.y);
    [TheGrid.x,TheGrid.y]=AdcircCppForward(TheGrid.x,TheGrid.y,lo0,la0);
    TheGrid=belint(TheGrid);
    TheGrid=el_areas(TheGrid);
    % also assume we need to convert the initial drifter locations
    [xi,yi]=AdcircCppForward(xi,yi,lo0,la0);
end
%TheGrid.strtree=ComputeStrTree(TheGrid);

if isempty(dt)
    dt=duration(time(2)-time(1));
end
if isempty(idt)
    idt=1;
else
    if idt<1 || idt > length(time)
        error('idt not between 1 and length of velocity field');
    end
end
% check timestep
if dt<eps
    error('Timestep<eps???  You''re kidding, right???')
end

if isempty(t1d),t1d=time(1);end
if isempty(t2d),t2d=time(end);end

if t1d<time(1) || t1d > time(end)
    error('Input starting time outside of velocity field times.')
end
if t2d<time(1) || t2d > time(end)
    error('Input ending time outside of velocity field times.')
end

disp(['Tracking Start: '  datestr(t1d,0) ])
disp(['Tracking End  : '  datestr(t2d,0) ])

% Check sizes of input drogue positions
if isempty(xi) || isempty(yi)
    idx=randi(length(TheGrid.x),10,1);
    xi=TheGrid.x(idx);
    yi=TheGrid.y(idx);
else
    if ~all(size(xi)==size(yi))
        error('Sizes of initial drog position arrays must be equal.')
    end
end

[mdrog,ndrog]=size(xi);
xi=xi(:);
yi=yi(:);

% Attach element finding arrays to fem_grid_struct
if ~isfield(TheGrid,'A')  || ... 
   ~isfield(TheGrid,'B')  || ...
   ~isfield(TheGrid,'A0') || ...
   ~isfield(TheGrid,'T')
    TheGrid=belint(TheGrid);
    disp('   BELINT info added to fem_grid_struct')
end
if ~isfield(TheGrid,'ar')
    disp('   EL_AREAS info added to fem_grid_struct')
    TheGrid=el_areas(TheGrid);
end

% Locate initial positions in grid
% j will be the array that keeps track of the
% current element number of each drog.  A NaN will
% be used to indicate drog in-activity, either because
% the drog is initially out-of-bounds, or because the
% drogue has left the domain during tracking.
j=findelem(TheGrid,[xi yi]);

% Allocate space for time history of positions
tt=t1d:dt:t2d;
tt=tt(1:idt:length(tt));
xx=NaN*(ones(size(tt))'*ones(size(xi')))';
yy=xx;
uu=xx;
vv=xx;
dt=dt*idt;
Nt=length(tt);

% global fdiag
% fdiag='linear';
%fdiag='real';

% get initial conditions, forcing
xx(:,1)=xi;
yy(:,1)=yi;
[ut,vt,tidx]=vel_interp(TheGrid,xi,yi,j,V,time,t1d);
uu(:,1)=ut;
vv(:,1)=vt;

%Draw initial positions on screen
if isfield(options,'draw') && options.draw
    %axx=[ 3.7319e+05   3.7635e+05   4.5794e+06   4.5819e+06];
    plotbnd(TheGridll,'LineWidth',2)
    lcontour(TheGridll,'z',0,'Color','k');
    [xi2,yi2]=AdcircCppInverse(xi,yi,lo0,la0);
    hdrog_initpos=line(xi2,yi2,'LineStyle','none','Marker','o', ...
        'MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
    axis('equal')
    axis([-70.565      -70.436       41.356       41.431])
    grid on
    grid minor
    drawnow
    CLim([-1.5 1.5])
    colormap(jet(30))
end

% Loop over time;
tnow=t1d;
iter=1;
fprintf('Starting: %d  %s\n',iter,datestr(tnow,0))
xnew=xi;
ynew=yi;

NLive=length(xi);
hdrog1=[];
hdrog2=[];
hdrog3=[];
hdrog4=[];
hdrogdead=[];
hdrogstuck=[];
hc=[];

while tnow<t2d
    if isfield(options,'verbose')
        fprintf('   (%4d/%4d) Integrating from %s to %s\n',iter,NLive,datestr(tnow,0),datestr(tnow+dt,0))
    end
    
    xnow=xnew;
    ynow=ynew;
    
    igood=find(~isnan(j));
    NLive=length(igood);
    % If j contains NaN's, these drogues are have exited the
    % domain at some previous timestep.
    
    if isempty(igood)
        fprintf('All drogues eliminated.\n')
        break
    end
    
    % Extract drogues currently in domain
    jgood=j(igood);
    xgood=xnow(igood);
    ygood=ynow(igood);
    
    [xnext,ynext,jnext]=track4(TheGrid,jgood,xgood,ygood,V,time,tnow,dt);
    j(igood)=jnext;
    xnew(igood)=xnext;
    ynew(igood)=ynext;
    
    tnow=tnow+dt;
    if tnow >= t2d, break, end
    
    iter=iter+1;
    
    [ut,vt,tidx]=vel_interp(TheGrid,xnext,ynext,jnext,V,time,tnow);
    xx(:,iter)      =xnew;
    yy(:,iter)      =ynew;
    uu(igood,iter)  =ut;
    vv(igood,iter)  =vt;
    tt(iter)        =tnow;
    
    %Update positions on screen
    if options.draw
        %      delete([hdrog1(:); hdrog2(:); hdrog3(:); hdrog4(:); hdrog_initpos(:); hdrogdead(:); hdrogstuck(:)])
        delete([hdrog1(:); hdrog2(:); hdrog3(:); hdrog4(:); hdrogdead(:); hdrogstuck(:)])
        delete(hc)
        i0=max(1,iter-options.lag);
        [xxc1,yyc1]=AdcircCppInverse(xx(:,i0:iter)',yy(:,i0:iter)',lo0,la0);
        %hdrog1=line(xxc1,yyc1,'LineStyle','-', 'Marker','.',...
        %    'LineWidth',.25,'Color','b');
        hdrog1=line(xxc1,yyc1,'LineStyle','none', 'Marker','.','LineWidth',.25,'Color','b');        
        %      hdrog1=line(xx',yy','LineStyle','-','Marker','.', ...
        %                'LineWidth',.25,'Color','b');
        
        %       if iter>options.lag
        %           hdrog4=line(xx(:,iter-20:iter)',yy(:,iter-20:iter)','LineStyle','-', ...
        %                 'LineWidth',1,'Color','b');
        %           hdrog3=line(xx(:,iter-40:iter-20)',yy(:,iter-40:iter-20)','LineStyle','-', ...
        %                 'LineWidth',.5,'Color','b');
        %       end
        [xxc2,yyc2]=AdcircCppInverse(xnew,ynew,lo0,la0);
        hdrog2=line(xxc2,yyc2,'LineStyle','none','Marker','o', ...
            'MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');
        istuck=xnew==xx(:,iter-1) & ynew==yy(:,iter-1);
        hdrogstuck=line(xxc2(istuck),yyc2(istuck),'LineStyle','none','Marker','o', ...
            'MarkerSize',4,'MarkerFaceColor','y','MarkerEdgeColor','k');
        hdrogdead=line(xxc2(isnan(j)),yyc2(isnan(j)),'LineStyle','none','Marker','o', ...
            'MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
        zz=Z{'zeta'}(tidx,:);
        hc=colormesh2d(TheGridll,zz);
        drawnow
    end
    
end

% Prepare for return
fprintf('Ending: [%d %s]\n',iter,datestr(tnow,0))
if ConvToXY
    % invert the forward projection
    [xx,yy]=AdcircCppInverse(xx,yy,lo0,la0);
end
D.tt=tt;
if mdrog==1
    D.xx=reshape(xx,ndrog,Nt);
    D.yy=reshape(yy,ndrog,Nt);
    D.uu=reshape(uu,ndrog,Nt);
    D.vv=reshape(vv,ndrog,Nt);
elseif ndrog==1
    D.xx=reshape(xx,mdrog,Nt);
    D.yy=reshape(yy,mdrog,Nt);
    D.uu=reshape(uu,mdrog,Nt);
    D.vv=reshape(vv,mdrog,Nt);
else
    D.xx=reshape(xx,mdrog,ndrog,Nt);
    D.yy=reshape(yy,mdrog,ndrog,Nt);
    D.uu=reshape(uu,mdrog,ndrog,Nt);
    D.vv=reshape(vv,mdrog,ndrog,Nt);
end

end



% PRIVATE FUNCTIONS

function jnew=locate_drog(TheGrid,x,y,j)
    jnew=j;
    % See if drogs are still in previously known element
    idx=belel(TheGrid,j,[x y]);
    inotfound=find(idx==0);
    % Get new elements, if not in previously known element
    if ~isempty(inotfound)
        idx=inotfound;
        jnew(idx)=findelem(TheGrid,[x(idx) y(idx)]);
        %   jnew(idx)=FindElementsInStrTree(TheGrid,x(idx), y(idx));
    end
end

function [xnew,ynew,jnew]=track2(TheGrid,j,x,y,V,timevec,t,dt)
    dts=seconds(dt);

    % k1
    [uk1,vk1]=vel_interp(TheGrid,x,y,j,V,timevec,t);  % eval field at current location

    % k2
    xtemp=x+.5*dts*uk1;
    ytemp=y+.5*dts*vk1;
    jtemp=locate_drog(TheGrid,xtemp,ytemp,j);  % relocate in elements
    [uk2,vk2,tidx]=vel_interp(TheGrid,xtemp,ytemp,jtemp,V,timevec,t+.5*dt);

    xnew=x+dts*(uk1 + uk2)/2;
    ynew=y+dts*(vk1 + vk2)/2;

    % If NaN is in j, then drog has left domain.  
    % insert last known location into arrays
    jnew=locate_drog(TheGrid,xnew,ynew,j);
    inan=find(isnan(jnew));
    if ~isempty(inan)
       xnew(inan)=x(inan);
       ynew(inan)=y(inan);
    end
     
end

function [xnew,ynew,jnew]=track4(TheGrid,j,x,y,V,timevec,t,dt)
    dts=seconds(dt);

    % k1
    [uk1,vk1,tidx]=vel_interp(TheGrid,x,y,j,V,timevec,t);

    % k2
    xtemp=x+.5*uk1*dts;
    ytemp=y+.5*vk1*dts;
    jtemp=locate_drog(TheGrid,xtemp,ytemp,j);
    [uk2,vk2,tidx]=vel_interp(TheGrid,xtemp,ytemp,jtemp,V,timevec,t+.5*dt);

    % k3
    xtemp=x+.5*uk2*dts;
    ytemp=y+.5*vk2*dts;
    jtemp=locate_drog(TheGrid,xtemp,ytemp,jtemp);
    [uk3,vk3,tidx]=vel_interp(TheGrid,xtemp,ytemp,jtemp,V,timevec,t+.5*dt);

    % k4
    xtemp=x+uk3*dts;
    ytemp=y+vk3*dts;
    jtemp=locate_drog(TheGrid,xtemp,ytemp,jtemp);
    [uk4,vk4,tidx]=vel_interp(TheGrid,xtemp,ytemp,jtemp,V,timevec,t+dt);

    xnew=x+(uk1 + 2*uk2 + 2*uk3 + uk4)*dts/6;
    ynew=y+(vk1 + 2*vk2 + 2*vk3 + vk4)*dts/6;

    % If NaN is in j, then drog has left domain.  
    % insert last known location into arrays
    jnew=locate_drog(TheGrid,xnew,ynew,j);
    inan=find(isnan(jnew));
    if ~isempty(inan)
       xnew(inan)=x(inan);
       ynew(inan)=y(inan);
    end
end
    
function retval=belel(TheGrid,j,xylist)
    %BELEL - determine if points are in elements
    % BELEL
    tol=eps*10000000;
    phi=basis2d(TheGrid,xylist,j);
    test=phi>=-tol & phi<=1+tol;
    retval=all(test'==1);
end

function [u,v,tidx]=vel_interp(TheGrid,x,y,j,V,time,t)
    global dim3d

    % Get the velocities at this time
    % Temporal interpolation of velocity slices to this time.
    it1=find(t<=time);
    it1=it1(1);
    it2=find(t>=time);
    it2=it2(length(it2));
    if it1==it2
        tfac1=1;tfac2=0;
    else
        tfac1=(time(it1)-t)/(time(it1)-time(it2));
        tfac2=1-tfac1;
    end
    if dim3d>0
        u1=V{'u-vel3D'}(it1,dim3d,:)';
        u2=V{'u-vel3D'}(it2,dim3d,:)';
        v1=V{'v-vel3D'}(it1,dim3d,:)';
        v2=V{'v-vel3D'}(it2,dim3d,:)';
    else
        % JML: add perturbations to u,v here
        u1=V{'u-vel'}(it1,:)';
        u2=V{'u-vel'}(it2,:)';
        v1=V{'v-vel'}(it1,:)';
        v2=V{'v-vel'}(it2,:)';
    end
    % diagnostic flows for debugging integrator
    % global fdiag
    % switch fdiag
    %     case 'linear'
    %         u=1.*ones(size(x));
    %         v=0.*ones(size(x));
    %     case 'circular'
    %         omega=2*pi/43200;
    %         a=.01;
    %         tt=seconds(t-time(1));
    %         u=-a.*cos(omega*tt)*ones(size(x));
    %         v= a.*sin(omega*tt)*ones(size(x));
    %     otherwise  %  the real flow
    % Depending on the number of particles to track, the
    % interpolation to (x,y,t) is done one of two ways.
    % It is not obvious, but the flop savings can be huge.
    if length(j)>150
        % Interpolate in time first,...
        uu=tfac1*u2 + tfac2*u1;
        vv=tfac1*v2 + tfac2*v1;
        % ... then, space
        u=interp_scalar(TheGrid,uu,x,y,j);
        v=interp_scalar(TheGrid,vv,x,y,j);
    else
        % Interpolate in space first, at the 2 time levels
        uu1=interp_scalar(TheGrid,u1,x,y,j);
        vv1=interp_scalar(TheGrid,v1,x,y,j);
        uu2=interp_scalar(TheGrid,u2,x,y,j);
        vv2=interp_scalar(TheGrid,v2,x,y,j);
        % Then, interpolate BETWEEN time levels
        u=tfac1*uu2 + tfac2*uu1;
        v=tfac1*vv2 + tfac2*vv1;
    end
    tidx=it1;
    % end
 end

%
%        Brian O. Blanton
%        RENCI
%        University of North Carolina
%        Chapel Hill, NC
%
%        brian_blanton@renci.org
%
%        March 2018
%
