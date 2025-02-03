%.........................................................................
%
%  Script to compute momentum balances from ADCIRC output in streamwise,
%  normal (s,n) coordinates
%
%  The sign convention used follows the right hand rule, i.e., the positive 
%  n direction is 90 deg to left of s.  Centrifugal accleration is positive
%  when curvature is concave to the left of the flow direction and negative
%  when curvature is concave to the right of the flow direction.
%
%  -g dzeta/dn and -g dzeta/ds are assumed to be on the RHS of their 
%  corresponding momentum equations.  All other terms are moved to the LHS 
%  of their momentum equations to emphasize the fact that these terms add 
%  together to produce the elevation gradients
% 
% Requires nctoolbox and adcirc_util libraries
%
%                     2/8/2022   - Rick Luettich
%                                   Brian Blanton
%                                   Matt Bilskie
%                    9/28/2022   - Matt Bilskie found a bug in CPP projection
% 
% 30 Jan 2022, orig, RL
% 01 Feb 2022, added url option, some error checking, BB
% 04 Feb 2022, added spherical coord corrections
% 05 Feb 2022, cleaned up code, corrected bugs in Bernoulli accel & dalpha/dt
% 06 Feb 2022, added time step inc as input, modified plots to work for 
%                                                    multiple time steps  
% 08 Feb 2022, added reading nodal attritsr_numbutes from fort.13, and inclusion of  
%                    spatially varying Mannings n in bottom stress calc.
% 09 Sept 2022, bug fix by Matt Bilskie on conversion from lon,lat to m?
% 18 Dec 2022, Matt's bug fix implemented and the sign of the terms
%                    adjusted as described above; plot surface elevation
% 19 Dec 2022, added water velocity vectors and wind velocity vectors; 
%                    revised multi-panel figure to 4x4
%

function momentum_sn(ts0,ts1,tsinc,url,varargin)
%function [g, times, dusdt, bernaccel, gdzetads, botfric, gdatmospresds,...
%    radstressgrads, LatStrgrads, totalmombals, ...
%    usdalphadt, centaccel, coraccel, gdzetadn, gdatmospresdn,...
%    windstressn, windstresss radstressgradn, LatStrgradn, totalmombaln  ...
%    twodprimcont] = momentum_sn(ts0,ts1,tsinc,url)

echo momentum_sn off

% ts0 - first timestep in ADCIRC output array to analyze
% Note to compute the time derivatives (centered aross current time step), 
%   it is necessary to start with times step 2
%
% For H. Michael NOPP simulations
%   COAMPS-TC
%        ts0=450;  = 15:00 hrs on 10/10
%        ts0=467;  = 17:50 hrs on 10/10 ~ max surge
%        output every 10 min
%  GAHM+NAM
%         ts0=279  =  15:00 hrs on 10/10
%         ts1=288  =  18:00 hrs on 10/10
%         output every 20 min
%
% ts1 - last timestep in ADCIRC output array to analyze
% Note to compute the time derivatives (centered aross current time step), 
%    it is necessary to end with the next to last time step
%
% tsinc - increment between timesteps to plot (i.e., tsinc=1 plot every time
%         step, tsinc=2 plot every other time step). 
%
% url - allows input files to be accessed across the network, e.g., from a
%       THREDDS server.  If using current MATLAB directory use '.'
%
% vargin - panel or animation
%
% needed constants.  Note phi0, lamda0, Manningn_def, Cfric_base,
% Cdrag_cap, ELSM, SmagMax, SmagMin are run dependent and should be modified
% to match values in the ADCIRC input file

grav=9.81;
rhow=1000.;          %kg/m^3
rhoa=1.15;           %kg/m^3
rhoa_o_rhow=rhoa/rhow;
atm_background=(1013*100)/(grav*rhow);       %background atm pressure in m H2O
REarth=6378206.4;    %Radius of the Earth (m) used in CPP spherical coord transformation
                     % hard coded into CPP conversion functions
omega = 7.29212e-5;  % Rate of earth rotation (1/s)
phi0=29.0;           %reference latitude for CPP spherical coord transformation (grid dependent)
                     % hard coded into ExtractGrid1
lamda0=265.5;        %reference longitude for CPP spherical coord transformation (grid dependent)
                     % hard coded into ExtractGrid1
Manningn_def=0.022;
Cfric_base=0.001;   % minimum bottom friction coeff
Cdrag_cap=0.0025;   % upper limit on wind friction coeff
ELSM=0.2;           % leading coefficient in Smagorinsky Lateral Viscosity
SmagMax=100;        % maximum allowable SmagEv value
SmagMin=0;          % minimum allowable SmagEv value

% process input file location
if nargin==0
    ts0=2;
    ts1=NaN;
    tsinc=1;
    url='./';
end
%if ~exit('tsinc')
%    tsinc=1;
%end
if ~exist('url')
    url='./';
end

panel=false;
animation=false;

k=1;
while k<length(varargin)+1
  switch lower(varargin{k})
    case 'panel'   
      panel=true;
      varargin([k])=[];    
    case 'animation'
      animation=true;
      varargin([k])=[];
    otherwise
      k=k+1;
  end
end

% use nctoolbox to read in ADCIRC netCDF files
fprintf('Opening ADCIRC netCDF files:  \n')
files={[url '/fort.63.nc']
       [url '/fort.64.nc']
       [url '/fort.73.nc']
       [url '/fort.74.nc']
       [url '/rads.64.nc']};
nnames={'nc63','nc64','nc73','nc74','ncrads64'};
for i=1:length(nnames)
    fprintf('   %s ... \n',files{i})
    try 
        nc=ncgeodataset(files{i});
        eval(sprintf('%s=%s',nnames{i},'nc;'))
    catch 
        error('Failed to open %s.\n',files{i})
    end
end

nodalFile = 'fort.13';
try
    nf = read_adcirc_fort13(nodalFile);
catch
    error('Failed to open %s.\n',nodalFile)
end
    
fprintf(' Finished reading ADCIRC files.\n')

% extract time and convert to date & time format
fprintf('Process times... ')
times=nctime(nc63);
delts=seconds(times(2)-times(1));
times_str=datestr(times);      %convert to string to use on figures
fprintf(' Done...\n')
fprintf('Delta_t = %0.1f\n',delts)

% error checking on the ending and beginning times, compute the number of
% times to process
if isnan(ts1)
    ts1=length(times);
elseif ts1>length(times)
    error('Input ts1 exceeds length of time vector.')
end
if ts0<2 || ts0>length(times) || ts0 > ts1 
    error('Input ts0 is < 2 or exceeds length of time.')
end
tsrange=1+(ts1-ts0)/tsinc;             %total number time steps to analyze
tsrange=tsrange-mod(tsrange,1);

%use Brian's adcirc_util routines to extract grid related information
%type g to see all of the variables it contains
%g.nn = nn - number of nodes
%g.ne = ne - number of elements
%g.e = (ne,3) array - element table, nodes ordered counter clockwise around element 
%g.x = (nn) vector - nodal lamda coordinate (deg)
%g.y = (nn) vector - nodal phi coordinate (deg)
%g.z = (nn) vector - nodal depth (m) positive down
%g.dx = (ne,3) array - x-coord diff of 2 opposite nodes on ele
%                                                  [2-3,3-1,1-2] - not used
%g.dy = (ne,3) array - y-coord diff of 2 opposite nodes on ele
%                                                  [2-3,3-1,1-2] - not used
%g.A = (ne,3) array - x-coord diff of 2 opposite nodes on ele [3-2,1-3,2-1]
%                            - matches a1,a2,a3 in ADCIRC theory report
%g.B = (ne,3) array - y-coord diff of 2 opposite nodes on ele [2-3,3-1,1-2]
%                            - matches b1,b2,b3 in ADCIRC theory report
%g.ar = (ne) vector - element areas (deg^2) 
%These are all replicated in cartesian coordinates (m) with coresponding
%                                  variables with _cart appended

fprintf('Extracting grid info... ')
g=ExtractGrid1(nc63);
fprintf(' Done...\n')

% Multiply CPP Correction factor cos(phi0)/cos(phi) to B_cart terms to
% correct x-derivatives
e=g.e;
A_cart=g.A_cart;
B_cart=g.B_cart.*cos(deg2rad(phi0))./cos(deg2rad(g.y(e)));
% B_cart=g.B_cart.*cos(deg2rad(g.y(e)))/cos(deg2rad(phi0)); % error pointed out by Matt B.

% Set up a nodal table that lists the elements connected to each node.
% this is done using ComputeNodeElementAdjacency icne(:,maxelenei+1:end)=[];

% icne - matrix of size [NN x maxconnections] where maxconnections is
%                 the maximum number of elements that contain any node.
%       The elements that connect to a node
% c    - number of elements associated with each node ([NN x 1])

fprintf('Computing node/element table for grid ... ')
[icne, c]=ComputeNodeElementAdjacency(g);
maxelenei=max(c);

% Replace "0" entries in the node table with a final ficticious element # (g.ne+1)
icne(icne==0)=g.ne+1;
fprintf(' Done...\n')

% setup companion nodal arrays to icne that have A & B values instead of element #
fprintf('Setting up companion nodal arrays ... ')
A_cart_icne=zeros(g.nn,maxelenei);
B_cart_icne=A_cart_icne; % zeros(g.nn,maxelenei);
for i=1:g.nn
    for j=1:maxelenei
       enum=icne(i,j);
       if(enum<g.ne+1) 
          k=find(g.e(enum,:)==i);
          A_cart_icne(i,j)=A_cart(j,k);
          B_cart_icne(i,j)=B_cart(j,k);
       end
    end
end
fprintf(' Done...\n')

% compute Coriolis parameters at nodes
Corf_n = 2*sin(deg2rad(g.y))*omega;

% create local array of the elemental areas for later use
ar_cart_e=g.ar_cart;
ar_cart_e(g.ne+1)=0;

% create pointer arrays to water surface elevation and u,v depth averaged
% velocities, atmospheric pressure, wind velocity and radiational stress
zeta=nc63{'zeta'};
ubar=nc64{'u-vel'};
vbar=nc64{'v-vel'};
apre=nc73{'pressure'};    % atmos pressure (m of H2O)
winx=nc74{'windx'};
winy=nc74{'windy'};
radstressgrad_x=ncrads64{'radstress_x'};
radstressgrad_y=ncrads64{'radstress_y'};

% loop through the timesteps
fprintf('Preallocating a boat-load of big arrays ...')
dusdt=zeros(tsrange,g.nn);  
usdalphadt=zeros(tsrange,g.nn);
centaccel=zeros(tsrange,g.nn);  
bernaccel=zeros(tsrange,g.nn);
us=zeros(tsrange,g.nn);
coraccel=zeros(tsrange,g.nn);
botfric=zeros(tsrange,g.nn);
gdzetads=zeros(tsrange,g.nn);
gdzetadn=zeros(tsrange,g.nn);
windstresss=zeros(tsrange,g.nn);
windstressn=zeros(tsrange,g.nn);
gdatmospresds=zeros(tsrange,g.nn);
gdatmospresdn=zeros(tsrange,g.nn);
radstressgrads=zeros(tsrange,g.nn);
radstressgradn=zeros(tsrange,g.nn);
LatStrgrads=zeros(tsrange,g.nn);
LatStrgradn=zeros(tsrange,g.nn);
totalmombals=zeros(tsrange,g.nn);
totalmombaln=zeros(tsrange,g.nn);
twodprimcont=zeros(tsrange,g.nn);
Rs=zeros(tsrange,g.nn);
alpha=zeros(tsrange,g.nn);
centovercor=zeros(tsrange,g.nn);
SmagEv1_n_t=zeros(tsrange,g.nn);
SmagEv2_n_t=zeros(tsrange,g.nn);
% manningN = zeros(1,g.nn);

fprintf(' Done...\n')

% Fill in the Manning's n values
idx = find(contains(nf.att_names(:),'manning'));
manningN = nf.atts{idx};

%Loop over desired timestep range to compute momentum terms

fprintf('Computing momentum terms ...\n')
for i=ts0:tsinc:ts1
   ii=(i-ts0)/tsinc+1;
   if rem(i,1)==0
       fprintf('Time step # %5d (%5d)\n',ii,i);
   end

   %create nodal and elemental arrays of ADCIRC variables at this timestep 
   %need to do this to get variables out of pointer status
   zeta_n=zeta(i,:)';       %zeta value
   zeta_e=zeta_n(g.e);      %put into an array of ele x 3nodes
   zetatm1_n=zeta(i-1,:)';  %zeta value at previous timestep
   zetatp1_n=zeta(i+1,:)';  %zeta value at next timestep
   H_n=zeta_n+g.z;          %total water depth
   H_e=H_n(g.e);            %put into an array of ele x 3nodes
   ubar_n=ubar(i,:)';       %depth average u velocity
   ubar_e=ubar_n(g.e);      %put into an array of ele x 3nodes
   vbar_n=vbar(i,:)';       %depth average v velocity
   vbar_e=vbar_n(g.e);      %put into an array of ele x 3nodes
   ubartm1_n=ubar(i-1,:)';  %depth average u velocity at previous timestep
   vbartm1_n=vbar(i-1,:)';  %depth average v velocity at previous timestep
   ubartp1_n=ubar(i+1,:)';  %depth average u velocity at next timestep
   vbartp1_n=vbar(i+1,:)';  %depth average v velocity at next timestep
   apre_n=apre(i,:);        %atmospheric pressure (m of H20)
   apre_e=apre_n(g.e);      %put into an array of ele x 3nodes
   winx_n=winx(i,:)';       %x component of wind velocity
   winy_n=winy(i,:)';       %y component of wind velocity
   radstressgradx_n=radstressgrad_x(i,:)'./H_n;       %x component of wave radiation stress gradient
   radstressgrady_n=radstressgrad_y(i,:)'./H_n;       %y component of wave radiation stress gradient
%    Cfric_n=grav*Manningn_def^2./nthroot(H_n,3);       %compute the quadratic bottom friction coeff
   nroot = nthroot(H_n,3);
   Cfric_n=grav*manningN(:,1).^2./nroot;               %compute the quadratic bottom friction coeff
   Cfric_n=max(Cfric_n,Cfric_base);                   %apply specified lower limit
   wins_n=sqrt(winx_n.*winx_n+winy_n.*winy_n);        %calculate wind speed
   Cdrag_n=0.001*(0.75+0.067*wins_n);                 %calculate wind drag coeff using the Garratt drag law
   Cdrag_n=min(Cdrag_n,Cdrag_cap);                    %apply specified upper limit 
   windstressx_n= rhoa_o_rhow*Cdrag_n.*((wins_n.*winx_n)./H_n);  %compute x-direction wind stress 
   windstressy_n= rhoa_o_rhow*Cdrag_n.*((wins_n.*winy_n)./H_n);  %compute y-direction wind stress 
   % compute elemental values of cartesian spatial gradient times element 
   % area of zeta, H, ubar, vbar, apre
   ardzetadx_cart=(zeta_e(:,1).*B_cart(:,1)+zeta_e(:,2).*B_cart(:,2)+zeta_e(:,3).*B_cart(:,3))/2;
   ardzetady_cart=(zeta_e(:,1).*A_cart(:,1)+zeta_e(:,2).*A_cart(:,2)+zeta_e(:,3).*A_cart(:,3))/2;
   ardHdx_cart=(H_e(:,1).*B_cart(:,1)+H_e(:,2).*B_cart(:,2)+H_e(:,3).*B_cart(:,3))/2;
   ardHdy_cart=(H_e(:,1).*A_cart(:,1)+H_e(:,2).*A_cart(:,2)+H_e(:,3).*A_cart(:,3))/2;
   ardubardx_cart=(ubar_e(:,1).*B_cart(:,1)+ubar_e(:,2).*B_cart(:,2)+ubar_e(:,3).*B_cart(:,3))/2;
   ardubardy_cart=(ubar_e(:,1).*A_cart(:,1)+ubar_e(:,2).*A_cart(:,2)+ubar_e(:,3).*A_cart(:,3))/2;
   ardvbardx_cart=(vbar_e(:,1).*B_cart(:,1)+vbar_e(:,2).*B_cart(:,2)+vbar_e(:,3).*B_cart(:,3))/2;
   ardvbardy_cart=(vbar_e(:,1).*A_cart(:,1)+vbar_e(:,2).*A_cart(:,2)+vbar_e(:,3).*A_cart(:,3))/2;
   ardapredx_cart=(apre_e(:,1).*B_cart(:,1)+apre_e(:,2).*B_cart(:,2)+apre_e(:,3).*B_cart(:,3))/2;
   ardapredy_cart=(apre_e(:,1).*A_cart(:,1)+apre_e(:,2).*A_cart(:,2)+apre_e(:,3).*A_cart(:,3))/2;

   % compute elemental values of SmagEv and the lateral stress  
   SmagEv1=ELSM*sqrt((ardubardx_cart-ardvbardy_cart).^2+(ardubardy_cart+ardvbardx_cart).^2);        %Smag Eddy viscosity in ADCIRC
   SmagEv2=ELSM*sqrt(ardubardx_cart.^2+ardvbardy_cart.^2+0.5*(ardubardy_cart+ardvbardx_cart).^2);   %Smag Eddy viscosity on wikipedia
   SmagEv1=min(SmagEv1,SmagMax);        %apply upper cap
   SmagEv2=min(SmagEv2,SmagMax);        %apply upper cap
   SmagEv1=max(SmagEv1,SmagMin);        %apply lower limit
   SmagEv2=max(SmagEv2,SmagMin);        %apply lower limit
   Tauxx=SmagEv1.*ardubardx_cart./ar_cart_e(1:end-1);
   Tauyy=SmagEv1.*ardvbardy_cart./ar_cart_e(1:end-1);
   Tauxy=SmagEv1.*ardubardy_cart./ar_cart_e(1:end-1);                %asymetric form
   Tauyx=SmagEv1.*ardvbardx_cart./ar_cart_e(1:end-1);                %asymetric form
   %Tauxy=SmagEv1.*((ardubardy_cart+ardvbardx_cart)./ar_cart_e(1:end-1));  %symmetric form
   %Tauyx=Tauxy;                                                      %symmetric form
   
   %add a ficticious element at the end of each elemental array with zero
   %values to account for padding in icne array.  Done above for ar_cart_e.
   ardzetadx_cart(g.ne+1)=0;
   ardzetady_cart(g.ne+1)=0;
   ardHdx_cart(g.ne+1)=0;
   ardHdy_cart(g.ne+1)=0;
   ardubardx_cart(g.ne+1)=0;   
   ardubardy_cart(g.ne+1)=0;
   ardvbardx_cart(g.ne+1)=0;   
   ardvbardy_cart(g.ne+1)=0;
   ardapredx_cart(g.ne+1)=0;
   ardapredy_cart(g.ne+1)=0;
   SmagEv1(g.ne+1)=0;
   SmagEv2(g.ne+1)=0;
   Tauxx(g.ne+1)=0;
   Tauyy(g.ne+1)=0;
   Tauxy(g.ne+1)=0;
   Tauyx(g.ne+1)=0;
         
   %aggregate derivatives at nodes - fast way
   
   dzetadx_cart_n=sum(ardzetadx_cart(icne),2);    %note the 2 means to sum across the columns (2nd array dimension)
   dzetady_cart_n=sum(ardzetady_cart(icne),2);
   dHdx_cart_n=sum(ardHdx_cart(icne),2);
   dHdy_cart_n=sum(ardHdy_cart(icne),2);
   dubardx_cart_n=sum(ardubardx_cart(icne),2);
   dubardy_cart_n=sum(ardubardy_cart(icne),2);   
   dvbardx_cart_n=sum(ardvbardx_cart(icne),2);
   dvbardy_cart_n=sum(ardvbardy_cart(icne),2);   
   dapredx_cart_n=sum(ardapredx_cart(icne),2);
   dapredy_cart_n=sum(ardapredy_cart(icne),2);
   SmagEv1_n=sum(ar_cart_e(icne).*SmagEv1(icne),2);
   SmagEv2_n=sum(ar_cart_e(icne).*SmagEv2(icne),2);
   LatStrxxdphidx_n=sum(Tauxx(icne).*B_cart_icne(icne),2)/2;  
   LatStryxdphidy_n=sum(Tauyx(icne).*A_cart_icne(icne),2)/2;
   LatStrxydphidx_n=sum(Tauxy(icne).*B_cart_icne(icne),2)/2;
   LatStryydphidy_n=sum(Tauyy(icne).*A_cart_icne(icne),2)/2;     
   ar_cart_n=sum(ar_cart_e(icne),2);
      
   % finish nodal derivative calc by dividing by the total elemental area
   % around each node
   dzetadx_cart_n=dzetadx_cart_n./ar_cart_n;
   dzetady_cart_n=dzetady_cart_n./ar_cart_n;
   dHdx_cart_n=dHdx_cart_n./ar_cart_n;
   dHdy_cart_n=dHdy_cart_n./ar_cart_n;
   dubardx_cart_n=dubardx_cart_n./ar_cart_n;
   dubardy_cart_n=dubardy_cart_n./ar_cart_n;
   dvbardx_cart_n=dvbardx_cart_n./ar_cart_n;
   dvbardy_cart_n=dvbardy_cart_n./ar_cart_n;
   dapredx_cart_n=dapredx_cart_n./ar_cart_n;
   dapredy_cart_n=dapredy_cart_n./ar_cart_n;
   SmagEv1_n=SmagEv1_n./ar_cart_n;
   SmagEv2_n=SmagEv2_n./ar_cart_n;
   LatStrxxdphidx_n=-3*LatStrxxdphidx_n./ar_cart_n;   %-3 comes from FE implementation and integration by parts
   LatStryxdphidy_n=-3*LatStryxdphidy_n./ar_cart_n;  
   LatStrxydphidx_n=-3*LatStrxydphidx_n./ar_cart_n;  
   LatStryydphidy_n=-3*LatStryydphidy_n./ar_cart_n;
   LatStrgradx_n=LatStrxxdphidx_n+LatStryxdphidy_n;
   LatStrgrady_n=LatStrxydphidx_n+LatStryydphidy_n;
   SmagEv1_n_t(ii,:)=SmagEv1_n;
   SmagEv2_n_t(ii,:)=SmagEv2_n;
   
   %compute stuff in x,y coordinates
   dudxpdvdy=dubardx_cart_n+dvbardy_cart_n;  
   
   %2D primitive continuity  dzeta/dt + d(UH)/dx + d(VH)/dy = 0
   %   dzeta/dt + UdH/dx + VdH/dy + H(dU/dx + dV/dy) = 0
   twodprimcont=zetatp1_n - zetatm1_n +(2*delts)*(ubar_n.*dHdx_cart_n + vbar_n.*dHdy_cart_n + H_n.*dudxpdvdy);
   
   %compute terms in s - n coordinates
   uss_n=ubar_n.*ubar_n+vbar_n.*vbar_n;    %along stream velocity^2
   us_mag_n=sqrt(uss_n);
   alpha(ii,:)=(atan2(ubar_n,vbar_n))';    %angle cw from N (radians)
   cosalpha_n=ubar_n./us_mag_n;
   sinalpha_n=vbar_n./us_mag_n;
   us(ii,:)=ubar_n.*cosalpha_n+vbar_n.*sinalpha_n;
 
   ustm1_mag_n=sqrt(ubartm1_n.*ubartm1_n+vbartm1_n.*vbartm1_n);
   cosalphatm1_n=ubartm1_n./ustm1_mag_n;
   sinalphatm1_n=vbartm1_n./ustm1_mag_n;   
   ustm1_n(:)=ubartm1_n.*cosalphatm1_n+vbartm1_n.*sinalphatm1_n;
   ustp1_mag_n=sqrt(ubartp1_n.*ubartp1_n+vbartp1_n.*vbartp1_n);
   cosalphatp1_n=ubartp1_n./ustp1_mag_n;
   sinalphatp1_n=vbartp1_n./ustp1_mag_n;    
   ustp1_n(:)=ubartp1_n.*cosalphatp1_n+vbartp1_n.*sinalphatp1_n;

   dusdt(ii,:)=(ustp1_n-ustm1_n)/(2*delts);       %centered difference
   alphatm1_n=(atan2(vbartm1_n,ubartm1_n))';      %angle ccw from E (radians)
   alphatp1_n=(atan2(vbartp1_n,ubartp1_n))';      %angle ccw from E (radians)
   dalpha_n=alphatp1_n-alphatm1_n;
   dalpha_n_index=dalpha_n>pi;
   dalpha_n=dalpha_n-2*pi*dalpha_n_index;
   dalpha_n_index=dalpha_n<-pi;
   dalpha_n=dalpha_n+2*pi*dalpha_n_index;
   usdalphadt(ii,:)=(us(ii,:).*dalpha_n)/(2*delts); %centered difference
   
   bernaccel(ii,:)=(ubar_n.*dubardx_cart_n+vbar_n.*dubardy_cart_n).*cosalpha_n + (ubar_n.*dvbardx_cart_n+vbar_n.*dvbardy_cart_n).*sinalpha_n;  %LHS of mom eq
   centaccel(ii,:)=(ubar_n.*dvbardx_cart_n+vbar_n.*dvbardy_cart_n).*cosalpha_n - (ubar_n.*dubardx_cart_n+vbar_n.*dubardy_cart_n).*sinalpha_n;  %LHS of mom eq
   coraccel(ii,:)=Corf_n.*us(ii,:)';                                                   %LHS of mom eq
   gdzetads(ii,:)=-grav*(dzetadx_cart_n.*cosalpha_n+dzetady_cart_n.*sinalpha_n);       %RHS of mom eq
   gdzetadn(ii,:)=-grav*(dzetady_cart_n.*cosalpha_n-dzetadx_cart_n.*sinalpha_n);       %RHS of mom eq
   gdatmospresds(ii,:)=grav*(dapredx_cart_n.*cosalpha_n+dapredy_cart_n.*sinalpha_n);   %LHS of mom eq
   gdatmospresdn(ii,:)=grav*(dapredy_cart_n.*cosalpha_n-dapredx_cart_n.*sinalpha_n);   %LHS of mom eq
   botfric(ii,:)=(Cfric_n.*uss_n)./H_n;                                                %LHS of mom eq
   windstresss(ii,:)=-(windstressx_n.*cosalpha_n+windstressy_n.*sinalpha_n);           %LHS of mom eq
   windstressn(ii,:)=-(windstressy_n.*cosalpha_n-windstressx_n.*sinalpha_n);           %LHS of mom eq
   radstressgrads(ii,:)=-(radstressgradx_n.*cosalpha_n+radstressgrady_n.*sinalpha_n);  %LHS of mom eq
   radstressgradn(ii,:)=-(radstressgrady_n.*cosalpha_n-radstressgradx_n.*sinalpha_n);  %LHS of mom eq
   LatStrgrads(ii,:)=-(LatStrgradx_n.*cosalpha_n+LatStrgrady_n.*sinalpha_n);           %LSH of mom eq
   LatStrgradn(ii,:)=-(LatStrgrady_n.*cosalpha_n-LatStrgradx_n.*sinalpha_n);           %LSH of mom eq
   Rs(ii,:)=(uss_n./centaccel(ii,:)')';
   centovercor(ii,:)=centaccel(ii,:)./coraccel(ii,:);
   totalmombals(ii,:)=dusdt(ii,:)+bernaccel(ii,:)-gdzetads(ii,:)+botfric(ii,:)+gdatmospresds(ii,:)+windstresss(ii,:)+radstressgrads(ii,:)+LatStrgrads(ii,:);
   totalmombaln(ii,:)=usdalphadt(ii,:)+centaccel(ii,:)+coraccel(ii,:)-gdzetadn(ii,:)+gdatmospresdn(ii,:)+windstressn(ii,:)+radstressgradn(ii,:)+LatStrgradn(ii,:);
   alpha(ii,:)=rad2deg(alpha(ii,:));
   zetap(ii,:)=zeta_n;
   ubarp(ii,:)=ubar_n;
   vbarp(ii,:)=vbar_n;
   usbarp(ii,:)=us_mag_n;
   winxp(ii,:)=winx_n;
   winyp(ii,:)=winy_n;
   winsp(ii,:)=wins_n;
   aprep(ii,:)=apre_n-atm_background;
end
fprintf(' Done.\n')

% set up for plots

%xmin=-85.9;       %H. Michael
%xmax=-84.7;       %H. Michael
%ymin=29.1;        %H. Michael
%ymax=30.3;        %H. Michael
%xmin=-85.0;       %H. Ian outer zoom
%xmax=-80.0;       %H. Ian outer zoom
%ymin=25.0;        %H. Ian outer zoom
%ymax=30.0;        %H. Ian outer zoom
xmin=-82.5;       %H. Ian inner zoom
xmax=-81.7;       %H. Ian inner zoom
ymin=26.1;        %H. Ian inner zoom
ymax=26.9;        %H. Ian inner zoom

AxisLims=[ xmin xmax ymin ymax];    %plot extents

% create a rectangular grids for vector plotting for wind and water velocity
% create a mask to zero out any values that lie outside the ADCIRC domain

fprintf('setting up rectangular grid points to plot velocity vectors...');

%delx=0.1;   %H. Michael
%dely=0.1;   %H. Michael
%delx=0.5;   %H. Ian outer zoom
%dely=0.5;   %H. Ian outer zoom
delx=0.08;   %H. Ian inner zoom
dely=0.08;   %H. Ian inner zoom

[xR,yR]=meshgrid(xmin:delx:xmax, ymin:dely:ymax);   % create regular grid to interpolate onto
XV=g.x(g.e);
YV=g.y(g.e);
IN = inpolygon(xR,yR,XV(1,:),YV(1,:));        %do once to initialize the IN matrix
for i=2:g.ne
   IN = IN + inpolygon(xR,yR,XV(i,:),YV(i,:));
end
fprintf(' Done.\n')

%   wmaxmag=100;       %H. Michael
%   wnumpoints=25;    %H. Michael
%   vmaxmag=4;        %H. Michael
%   vnumpoints=20;    %H. Michael
%   wmaxmag=75;       %H. Ian outer zoom
%   wnumpoints=50;    %H. Ian outer zoom
%   vmaxmag=2;        %H. Ian outer zoom
%   vnumpoints=50;    %H. Ian outer zoom
   wmaxmag=100;       %H. Ian inner zoom
   wnumpoints=50;     %H. Ian inner zoom
   vmaxmag=4;         %H. Ian inner zoom
   vnumpoints=20;     %H. Ian inner zoom

if panel

   titles1={
       'wind velocity (m/s)'
      'atmospheric pressure deficit (m H2O)'
       'water velocity (m/s)'
       'elevation (m) NAVD88'
       '-g dzeta/dn (m/s^2)'
       'us dalpha/dt (m/s^2)'
       'Centrifugal Acceleration (m/s^2)'
       'Coriolis Acceleration (m/s^2)'
       'wind stress/H n (m/s^2)'
       'g datmospres/dn (m/s^2)'
       'wave radiation stress grad/H n (m/s^2)'
       'lateral stress grad/H n (m/s^2)'
       'total n mom bal (m/s^2)'
       'abs(Centrifugal/Coriolis)'
       'Radius of curvature (km)'
       'alpha (deg cw from North)'
       };
   titles2={
       'wind velocity (m/s)'
       'atmospheric pressure deficit (m H2O)'
       'water velocity (m/s)'
       'elevation (m) NAVD88'
       '-g dzeta/ds (m/s^2)'
       'dus/dt (m/s^2)'
       'Bernoulli Acceleration (m/s^2)'
       'Bottom Stress/H (m/s^2)'
       'wind stress/H s (m/s^2)'
       'g datmospres/ds (m/s^2)'
       'wave radiation stress grad/H s (m/s^2)'
       'lateral stress grad/H s (m/s^2)'
       'total s mom bal (m/s^2) '
       'SmagEv1 (m^2/s)'
       'SmagEv2 (m^2/s)'
       'mass conservation error (m)'
       };

   plotcommands1={
      ['colormesh2d(g,winsp(tsa_num,:)''); hold; ' ...
       'curvvec(xR,yR,windxR,windyR,''color'',''k'',''minmag'',0,''maxmag'',wmaxmag,''numpoints'',wnumpoints,''thin'',1);']    
      ['colormesh2d(g,aprep(tsa_num,:)''); hold; ' ...
       'curvvec(xR,yR,windxR,windyR,''color'',''k'',''minmag'',0,''maxmag'',wmaxmag,''numpoints'',wnumpoints,''thin'',1);'] 
      ['colormesh2d(g,usbarp(tsa_num,:)''); hold; ' ...
       'curvvec(xR,yR,ubarR,vbarR,''color'',''k'',''minmag'',0,''maxmag'',vmaxmag,''numpoints'',vnumpoints,''thin'',1);']
       'colormesh2d(g,zetap(tsa_num,:)'');'
       'colormesh2d(g,gdzetadn(tsa_num,:)'');'
       'colormesh2d(g,usdalphadt(tsa_num,:)'');'
       'colormesh2d(g,centaccel(tsa_num,:)'');'
       'colormesh2d(g,coraccel(tsa_num,:)'');'
       'colormesh2d(g,windstressn(tsa_num,:)'');'
       'colormesh2d(g,gdatmospresdn(tsa_num,:)'');'
       'colormesh2d(g,radstressgradn(tsa_num,:)'');'
       'colormesh2d(g,LatStrgradn(tsa_num,:)'');'
       'colormesh2d(g,totalmombaln(tsa_num,:)'');'
       'colormesh2d(g,centovercor(tsa_num,:)'');'
       'colormesh2d(g,Rs(tsa_num,:)''/1000);'
       'colormesh2d(g,alpha(tsa_num,:)'');'
       };
   plotcommands2={
       ['colormesh2d(g,winsp(tsa_num,:)''); hold; ' ...
       'curvvec(xR,yR,windxR,windyR,''color'',''k'',''minmag'',0,''maxmag'',wmaxmag,''numpoints'',wnumpoints,''thin'',1);']
       ['colormesh2d(g,aprep(tsa_num,:)''); hold; ' ...
       'curvvec(xR,yR,windxR,windyR,''color'',''k'',''minmag'',0,''maxmag'',wmaxmag,''numpoints'',wnumpoints,''thin'',1);'] 
       ['colormesh2d(g,usbarp(tsa_num,:)''); hold; ' ...
       'curvvec(xR,yR,ubarR,vbarR,''color'',''k'',''minmag'',0,''maxmag'',vmaxmag,''numpoints'',vnumpoints,''thin'',1);']
       'colormesh2d(g,zetap(tsa_num,:)'');'
       'colormesh2d(g,gdzetads(tsa_num,:)'');'
       'colormesh2d(g,dusdt(tsa_num,:)'');'
       'colormesh2d(g,bernaccel(tsa_num,:));'
       'colormesh2d(g,botfric(tsa_num,:)'');'
       'colormesh2d(g,windstresss(tsa_num,:)'');'
       'colormesh2d(g,gdatmospresds(tsa_num,:)'');'
       'colormesh2d(g,radstressgrads(tsa_num,:)'');'
       'colormesh2d(g,LatStrgrads(tsa_num,:)'');'
       'colormesh2d(g,totalmombals(tsa_num,:)'');'
       'colormesh2d(g,SmagEv1_n_t(tsa_num,:)'');'
       'colormesh2d(g,SmagEv2_n_t(tsa_num,:)'');'
       'colormesh2d(g,twodprimcont(:)'');'
       };

   fignum=0;

   for it=ts0:tsinc:ts1
  
       tsa_num=(it-ts0)/tsinc+1;
    
       windxR=griddata(g.x,g.y,winxp(tsa_num,:),xR,yR).*IN;
       windyR=griddata(g.x,g.y,winyp(tsa_num,:),xR,yR).*IN;
       ubarR=griddata(g.x,g.y,ubarp(tsa_num,:),xR,yR).*IN;
       vbarR=griddata(g.x,g.y,vbarp(tsa_num,:),xR,yR).*IN;

       fprintf('Creating part 1 mosaic figure # %5d (%5d)\n',tsa_num,it);
       close('all')

       fignum=fignum+1;
       fig=figure(fignum);     
       tiledLayout=tiledlayout(4,4);

       set(gcf,'WindowState','maximized')
    
       for axnum=1:length(plotcommands1)
           nti=nexttile(axnum);        
           nti.XTickLabels={};
           nti.YTickLabels={};
        
           eval(plotcommands1{axnum})
           axis('equal')
           axis(AxisLims)
           colormap(jet(20))
           switch axnum
               case 1                       %wind speed (m/s)
%                   caxis([0 70])            %H. Michael     
                   caxis([0 60])            %H. Ian
               case 2                       %atm press deficit (m H2O)
%                   caxis([-1 0])            %H. Michael
                   caxis([-0.8 0])          %H. Ian
               case 3                       %water speed (m/s)
%                   caxis([0 3])             %H. Michael
                   caxis([0 2.5])           %H. Ian            
               case 4                       %elevation
%                   caxis([-1 4])            %H. Michael
                   caxis([-1 3])            %H. Ian
%               case 5                       %-g dzeta/dn
%               case 6                       % us dalpha/dt
%               case 7                       % centrifugal accel
%                   caxis([-0.0005 0.0005])
%               case 8                       % Coriolis accel
%                   caxis([-0.0005 0.0005])
%               case 9                       % wind stress n
%               case 10                      % atm press grad n
%               case 11                      % rad stress n
               case 12                      % lat stress n
                   caxis([-0.000001 0.000001])
               case 13                      %n-momentum balance
                   caxis([-0.0002 0.0002])
               case 14                      %ratio of centrifugal to coriolis accel
                   caxis([0 3])              
               case 15                      %radius of curviture
                   caxis([-100 100])
               case 16                       %alpha
                   caxis([-180 180]) 
               otherwise
                   caxis([-0.001 0.001])
           end

           title([int2str(axnum) ': ' titles1{axnum} ])
           sgtitle(datestr(times(it)),'FontSize',24) 
           colorbar('EastOutside')
       end
       saveas(fig,sprintf('sn_part1.%02d.png',it))
%      exportgraphics(fig, sprintf('sn.%02d.png',it), 'Resolution', 144);

       fignum=fignum+1;
       fig=figure(fignum);     
       tiledLayout=tiledlayout(4,4);
       fprintf('Creating part 2 mosaic figure # %5d (%5d)\n',tsa_num,it);
       set(gcf,'WindowState','maximized')
    
       for axnum=1:length(plotcommands2)
           nti=nexttile(axnum);        
           nti.XTickLabels={};
           nti.YTickLabels={};
        
           eval(plotcommands2{axnum})
           axis('equal')
           axis(AxisLims)
           colormap(jet(20))
           switch axnum
               case 1                       %wind speed (m/s)
%                   caxis([0 70])            %H. Michael     
                   caxis([0 60])            %H. Ian
               case 2                       %atm press deficit (m H2O)
%                   caxis([-1 0])            %H. Michael
                   caxis([-0.8 0])          %H. Ian
               case 3                       %water speed (m/s)
%                   caxis([0 3])             %H. Michael
                   caxis([0 2.5])           %H. Ian            
               case 4                       %elevation
%                   caxis([-1 4])            %H. Michael
                   caxis([-1 3])            %H. Ian
%               case 5                       %-g dzeta/dn
%               case 6                       % dus/dt
%               case 7                       % bernulli acceleration
%               case 8                       % bottom friction
%               case 9                       % wind stress s
%               case 10                      % atm press grad s
%               case 11                      % rad stress s
               case 12                      % lat stress s
                   caxis([-0.000001 0.000001])
               case 13                      % s-momentum balance
                   caxis([-0.0002 0.0002])
               case {14, 15}                 %Smag EV
                   caxis([0 100])
               case 16                       %mass conservation
                   caxis([-0.5 0.5])
               otherwise
                   caxis([-0.001 0.001])
           end

           title([int2str(axnum) ': ' titles2{axnum} ])
           sgtitle(datestr(times(it)),'FontSize',24) 
           colorbar('EastOutside')
       end
       saveas(fig,sprintf('sn_part2.%02d.png',it))
%      exportgraphics(fig, sprintf('sn.%02d.png',it), 'Resolution', 144);

   end
end                   %end of plot panel
    
% plot animations of individual terms
elevation=false;
velocity=false;
centrifugal=false;
corriolis=false;
centocorr=false;
dzetadn=false;
wstressn=false;
rstressn=false;
dpresdn=false;
mombaln=false;
dzetads=false;
bernoulli=false;
bstress=false;
wstresss=false;
rstresss=false;
dpresds=false;
mombals=false;

if animation
   elevation=true;
%   velocity=true;
   centrifugal=true;
   corriolis=true;
   centocorr=true;
   dzetadn=true;
   wstressn=true;
   rstressn=true;
   dpresdn=true;
   mombaln=true;
   dzetads=true;
   bernoulli=true;
   bstress=true;
   wstresss=true;
   rstresss=true;
   dpresds=true;
   mombals=true;

   fignum=99;
   for it=ts0:tsinc:ts1
   
      tsa_num=(it-ts0)/tsinc+1;
      tsr_num=it;   

      fprintf('Creating individual term figures # %5d (%5d)\n',tsa_num,it);
      close('all')
      nametime=datestr(times(tsr_num),'DDmmmYYYY hhMMSS_');
      
      if elevation

         fignum=fignum+1;
         fig=figure(fignum);  
         axis(AxisLims);            
         colormesh2d(g,zetap(tsa_num,:)');
         colormap(jet(20))
         caxis([-1 4])
         colorbar
         title(['elevation (m) NAVD88 ' times_str(tsr_num,:)]);
         nameroot = 'elevation';
%        saveas(fig,[nametime nameroot '.png'])
   
         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            e4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            e4.FrameRate=2;   %2 frames / second
            e4.Quality=100;
            open(e4);
         end
         writeVideo(e4,frames);
      end
   
      if velocity
         fignum=fignum+1;
         fig=figure(fignum);     %along stream velocity (m/s)
         axis(AxisLims);    %set plot extents

         ubarR=griddata(g.x,g.y,ubarp(tsa_num,:)',xR,yR).*IN;
         vbarR=griddata(g.x,g.y,vbarp(tsa_num,:)',xR,yR).*IN;
         ubarR(isnan(ubarR))=0;
         vbarR(isnan(vbarR))=0;
         
         colormesh2d(g,us(tsa_num,:)');
         hold;
         curvvec(xR,yR,ubarR,vbarR,'color','k','minmag',0,'maxmag',vmaxmag,'numpoints',vnumpoints,'thin',1);
         colormap(jet(20))
         clim([0 3])
         colorbar
         title(['us along stream vel (m/s) ' times_str(tsr_num,:)]);

         nameroot = 'velocity';         
%        saveas(fig,[nametime 'alongstreamvel.png'])

         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            v4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            v4.FrameRate=2;   %2 frames / second
            v4.Quality=100;
            open(v4);
         end
         writeVideo(v4,frames);
      end 

      if centrifugal

         fignum=fignum+1;
         fig=figure(fignum);     %centrifugal acceleration
         axis(AxisLims);    %set plot extents
         colormesh2d(g,centaccel(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.001 0.001])
         colorbar
         title(['Centrifugal Acceleration (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'cent_accel';
%        saveas(fig,[nametime 'n_centrifugalaccel.png'])
 
         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            cn4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            cn4.FrameRate=2;   %2 frames / second
            cn4.Quality=100;
            open(cn4);
         end
         writeVideo(cn4,frames);
      end

      if corriolis

         fignum=fignum+1;
         fig=figure(fignum);     %Coriolis acceleration
         axis(AxisLims);    %set plot extents
         colormesh2d(g,coraccel(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.001 0.001])
         colorbar
         title(['Coriolis Acceleration (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'corriolis_accel';
%        saveas(fig,[nametime 'n_coriolisaccel.png'])

         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            cr4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            cr4.FrameRate=2;   %2 frames / second
            cr4.Quality=100;
            open(cr4);
         end
         writeVideo(cr4,frames);
      end

      if centocorr    %centrifugal / coriolis  

         fignum=fignum+1;
         fig=figure(fignum);     
         axis(AxisLims);   
         colormesh2d(g,abs(centovercor(tsa_num,:)'));
         colormap(jet(20))
         caxis([0 3])
         colorbar
         title(['abs(Centrifugal/Coriolis) ' times_str(tsr_num,:)]);
         nameroot = 'centovercorr_accel';
%        saveas(fig,[nametime 'n_centovercoraccel.png'])
 
         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            coc4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            coc4.FrameRate=2;   %2 frames / second
            coc4.Quality=100;
            open(coc4);
         end
         writeVideo(coc4,frames);
      end

      if dzetadn

         fignum=fignum+1;
         fig=figure(fignum);    
         axis(AxisLims);    %set plot extents
         colormesh2d(g,gdzetadn(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.001 0.001])
         colorbar
         title(['-g dzeta/dn (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'gdzetadn';
%        saveas(fig,[nametime 'n_gdzetadn.png'])
 
         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            den4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            den4.FrameRate=2;   %2 frames / second
            den4.Quality=100;
            open(den4);
         end
         writeVideo(den4,frames);
      end

      if wstressn

         fignum=fignum+1;
         fig=figure(fignum);     %Surface winds stress/H in n direction
         axis(AxisLims);    %set plot extents
         colormesh2d(g,windstressn(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.001 0.001])
         colorbar
         title(['wind stress/H n (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'windstressn';
%        saveas(fig,[nametime 'n_windstressn.png'])

         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            wn4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            wn4.FrameRate=2;   %2 frames / second
            wn4.Quality=100;
            open(wn4);
         end
         writeVideo(wn4,frames);
      end

      if rstressn

         fignum=fignum+1;
         fig=figure(fignum);     %wave radiation stress grad/H in n direction
         axis(AxisLims);    %set plot extents
         colormesh2d(g,radstressgradn(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.001 0.001])
         colorbar
         title(['wave radiation stress/H n (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'radstressn';
%        saveas(fig,[nametime 'n_waveradstressgrad.png'])
  
         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            rn4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            rn4.FrameRate=2;   %2 frames / second
            rn4.Quality=100;
            open(rn4);
         end
         writeVideo(rn4,frames);
      end

      if dpresdn

         fignum=fignum+1;
         figure(fignum);     %g datmospres/dn
         axis(AxisLims);    %set plot extents
         colormesh2d(g,gdatmospresdn(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.001 0.001])
         colorbar
         title(['g datmospres/dn (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'dpresdn';
%         saveas(fig,[nametime 'n_gdatmospresdn.png'])

         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            dpn4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            dpn4.FrameRate=2;   %2 frames / second
            dpn4.Quality=100;
            open(dpn4);
         end
         writeVideo(dpn4,frames);
      end

      if mombaln
    
         fignum=fignum+1;
         fig=figure(fignum);  
         axis(AxisLims);    %set plot extents
         colormesh2d(g,totalmombaln(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.0002 0.0002])
         colorbar
         title(['total n mom bal (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'mombaln';
%       saveas(fig,[nametime 'n_totalmombal.png'])

         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            mbn4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            mbn4.FrameRate=2;   %2 frames / second
            mbn4.Quality=100;
            open(mbn4);
         end
         writeVideo(mbn4,frames);
      end

    if dzetads

         fignum=fignum+1;
         fig=figure(fignum);
         axis(AxisLims);    %set plot extents
         colormesh2d(g,gdzetads(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.001 0.001])
         colorbar
         title(['-g dzeta/ds (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'gdzetads';
%         saveas(fig,[nametime 's_gdzetads.png'])

         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            des4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            des4.FrameRate=2;   %2 frames / second
            des4.Quality=100;
            open(des4);
         end
         writeVideo(des4,frames);
      end 

      if bernoulli

         fignum=fignum+1;
         fig=figure(fignum);     %bernoulli acceleration
         axis(AxisLims);    %set plot extents
         colormesh2d(g,bernaccel(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.001 0.001])
         colorbar
         title(['Bernoulli Acceleration (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'bernoulli_accel';
%         saveas(fig,[nametime 's_bernoulliaccel.png'])

         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            ber4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            ber4.FrameRate=2;   %2 frames / second
            ber4.Quality=100;
            open(ber4);
         end
         writeVideo(ber4,frames);
      end 

      if bstress

         fignum=fignum+1;
         fig=figure(fignum);  
         axis(AxisLims);    %set plot extents
         colormesh2d(g,botfric(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.001 0.001])
         colorbar
         title(['Bottom Stress/H (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'botstress';
%        saveas(fig,[nametime 's_bottomstress.png'])

         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            bs4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            bs4.FrameRate=2;   %2 frames / second
            bs4.Quality=100;
            open(bs4);
         end
         writeVideo(bs4,frames);
      end

      if wstresss
 
         fignum=fignum+1;
         fig=figure(fignum);     %Surface winds stress/H in s direction
         axis(AxisLims);    %set plot extents
         colormesh2d(g,windstresss(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.001 0.001])
         colorbar
         title(['wind stress/H s (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'windstresss';
%         saveas(fig,[nametime 's_windstress.png'])

         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            ws4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            ws4.FrameRate=2;   %2 frames / second
            ws4.Quality=100;
            open(ws4);
         end
         writeVideo(ws4,frames);
      end

      if rstresss

         fignum=fignum+1;
         fig=figure(fignum);     %wave radiation stress grad/H in s direction
         axis(AxisLims);    %set plot extents
         colormesh2d(g,radstressgrads(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.001 0.001])
         colorbar
         title(['wave radiation stress grad/H s (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'radstresss';
%         saveas(fig,[nametime 's_waveradstressgrad.png'])

         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            rs4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            rs4.FrameRate=2;   %2 frames / second
            rs4.Quality=100;
            open(rs4);
         end
         writeVideo(rs4,frames);
      end 

      if dpresds
          
         fignum=fignum+1;
         fig=figure(fignum);     %g datmospres/ds
         axis(AxisLims);    %set plot extents
         colormesh2d(g,gdatmospresds(tsa_num,:)');
         colormap(jet(20))
         caxis([-0.001 0.001])
         colorbar
         title(['g datmospres/ds (m/s^2) ' times_str(tsr_num,:)]);
         nameroot = 'dpresds';
%         saveas(fig,[nametime 's_gdatmospresds.png'])

         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            dps4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            dps4.FrameRate=2;   %2 frames / second
            dps4.Quality=100;
            open(dps4);
         end
         writeVideo(dps4,frames);
      end

      if mombals

        fignum=fignum+1;
        fig=figure(fignum);     %total s momentum balance
        axis(AxisLims);    %set plot extents
        colormesh2d(g,totalmombals(tsa_num,:)');
        colormap(jet(20))
        caxis([-0.0002 0.0002])
        colorbar
        title(['total s mom bal (m/s^2)  ' times_str(tsr_num,:)]);
        nameroot = 'mombals';
%        saveas(fig,[nametime 'n_totalmombal.png'])

         frames = getframe(fig);

%        gif animation    
         im = frame2im(frames); 
         [A,maps] = rgb2ind(im,256);
         if tsa_num == 1
            imwrite(A,maps,[nameroot '.gif'],'gif','LoopCount',Inf,'DelayTime',0.5);
         else
            imwrite(A,maps,[nameroot '.gif'],'gif','WriteMode','append','DelayTime',0.5);
         end

%        MPEG-4 animation    
         if tsa_num == 1
            mbs4=VideoWriter([nameroot '.mp4'],'MPEG-4');
            mbs4.FrameRate=2;   %2 frames / second
            mb4.Quality=100;
            open(mbs4);
         end
         writeVideo(mbs4,frames);
      end



% fignum=fignum+1;
% fig=figure(fignum);    %d us / dt in s equation
% axis(AxisLims);    %set plot extents
% colormesh2d(g,dusdt(tsa_num,:)');
% colormap(jet(20))
% caxis([-0.001 0.001])
% colorbar
% title([' dus/dt (m/s^2) ' times_str(tsr_num,:)]);
% saveas(fig,[nametime 's_dusdt.png'])

% fignum=fignum+1;
% fig=figure(fignum);    %lateral stress grads/H in n direction
% axis(AxisLims);    %set plot extents
% colormesh2d(g,LatStrgradn(tsa_num,:)');
% colormap(jet(20))
% caxis([-0.001 0.001])
% colorbar
% title(['lateral stress grad/H n (m/s^2) ' times_str(tsr_num,:)]);
% saveas(fig,[nametime 'n_laterstressgrad.png'])

% fignum=fignum+1;
% fig=figure(fignum);     %Radius of curvature
% axis(AxisLims);    %set plot extents
% colormesh2d(g,Rs(tsa_num,:)'/1000);
% colormap(jet(20))
% caxis([-100 100])
% colorbar
% title(['Radius of curvature (km) ' times_str(tsr_num,:)]);
% saveas(fig,[nametime 'RadiusCurvature.png'])

% fignum=fignum+1;
% fig=figure(fignum);     %lateral stress gradient/H in s direction
% axis(AxisLims);    %set plot extents
% colormesh2d(g,LatStrgrads(tsa_num,:)');
% colormap(jet(20))
% caxis([-0.001 0.001])
% colorbar
% title(['lateral stress grad/H s (m/s^2) ' times_str(tsr_num,:)]);
% saveas(fig,[nametime 's_laterstressgrad.png'])

% fignum=fignum+1;
% fig=figure(fignum);     %along stream angle (deg cw from N)
% axis(AxisLims);    %set plot extents
% colormesh2d(g,alpha(tsa_num,:)');
% colormap(jet(20))
% caxis([-180 180])
% colorbar
% title(['alpha (deg cw from North) ' times_str(tsr_num,:)]);
% saveas(fig,[nametime 'alongstreamang.png'])
 
% fignum=fignum+1;
% fig=figure(fignum);     %SmagEv1 (m/s)
% axis(AxisLims);    %set plot extents
% colormesh2d(g,SmagEv1_n');
% colormap(jet(20))
% caxis([0 50])
% colorbar
% title(['SmagEv1 (m^2/s) ' times_str(tsr_num,:)]);
% saveas(fig,[nametime 'SmagEv1.png'])
 
% fignum=fignum+1;
% fig=figure(fignum);     %SmagEv2 (m/s)
% axis(AxisLims);    %set plot extents
% colormesh2d(g,SmagEv2_n');
% colormap(jet(20))
% caxis([0 50])
% colorbar
% title(['SmagEv2 (m^2/s) ' times_str(tsr_num,:)]);
% saveas(fig,[nametime 'SmagEv2.png'])

% fignum=fignum+1;
% fig=figure(fignum); %2D Primitive Continuity
% axis(AxisLims);    %set plot extents
% colormesh2d(g,twodprimcont(:)');
% colormap(jet(20))
% caxis([-0.5 0.5])
% colorbar
% title(['Mass Conservation Error (m) ' times_str(tsr_num,:)]);
% saveas(fig,[nametime 'MassConservation.png'])
 
% fignum=fignum+1;
% fig=figure(fignum); %dudx + dvdy
% axis(AxisLims);    %set plot extents
% colormesh2d(g,dudxpdvdy(:)');
% colormap(jet(20))
% caxis([-0.0001 0.0001])
% colorbar
% title(['du/dx + dv/dy (1/s) ' times_str(tsr_num,:)]);
% saveas(fig,[nametime 'dudxpdvdy.png'])

   end
end

