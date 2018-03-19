function elev=resynt_ha_elev(ap,w,t)
%RESYNT_HA_ELEV resynthesize harmonic analysis results from ADCIRC
% RESYNT_HA_ELEV recombines the harmonic analysis results from ADCIRC
% for at a vector of times.  
%
%  Inputs: ap - amp/phase data matrix as output by READ_ADCIRC_HA_ELEV
%          w  - vector of frequencies (rad/sec), one per amp/pha comp
%          t  - vector of times at which to compute elevation (seconds)
%
% Outputs: elev - array of elevation fields of size [nn X length(t)]
%                 where nn is length(ap) 
%
% Call as: elev=resynt_ha_elev(ap,freq,times);
% 
if nargin~=3
   error('RESYNT_HA_ELEV requires 3 input arguments')
end

[mt,nt]=size(t);
if nt*mt>mt &  nt*mt>nt
   error('Time vector to RESYNT_HA_ELEV must be 1-D')
end
[mw,nw]=size(w);
if nw*mw>mw &  nw*mw>nw
   error('Frequency vector to RESYNT_HA_ELEV must be 1-D')
end

w=w(:);
t=t(:);

[nn,nf]=size(ap);
nf=nf/2;

d2r=pi/180;
amp=ap(:,1:nf);
pha=ap(:,nf+1:nf*2)*d2r;

elev=NaN*ones(nn,length(t));

for i=1:length(t)
   phase_arg=w*ones(1,nn)*t(i)+pha';
   elev(:,i)=sum(amp'.*cos(phase_arg))';
end
