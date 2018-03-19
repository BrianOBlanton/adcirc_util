%OPNML ADCIRC ROUTINES
%OPNML ADCIRC ROUTINES
% GRD_TO_OPNML Convert an ADCIRC grd file to an OPNML fem_grid_struct.
%     Convert an ADCIRC grd file to an OPNML fem_grid_struct.
%     ADCIRC grid information assumed in "fort.14" format.
%     The boundary/island information at the tail of the fort.14
%     file is ignored.
% Call as: fem_grid_struct=grd_to_opnml(fort14name);
%
% READ_ADCIRC_HA_ELEV read ADCIRC elevation output file.
%     This routine reads in the contents of an ADCIRC elevation file,
%     as typically output from an harmonic analysis output to a fort.53
%     file.  fort.53 is the default output file for the global harmonic
%     analysis of the elevation field.  
% Call as: [DATA,FREQ,PERNAMES]=read_adcirc_ha_elev(fname,flag,gname);
%
% READ_ADCIRC_HA_VEL read ADCIRC velocity output file
%     This routine reads in the contents of an ADCIRC elevation file,
%     as typically output from an harmonic analysis output to a fort.54
%     file.  fort.54 is the default output file for the global harmonic
%     analysis of the velocity field.  
% Call as: [DATA,FREQ,PERNAMES]=read_adcirc_ha_vel(fname,flag,gname);
%
% READ_ADCIRC_FORT61 read ADCIRC station elevation output file
%  READ_ADCIRC_FORT61 reads in the contents of an ADCIRC station 
%  elevation file, typically output to IO unit 61 (fort.61).  
%  This output file is described in the ADCIRC Users Manual as 
%  "Elevation Time Series at Specified Elevation Recording Stations (fort.61)". 
%  Time is converted from ADCIRC time to Gregorian time.
%  Call as: D=read_adcirc_fort61(fname,asd);
%
% READ_ADCIRC_FORT62 read ADCIRC station velocity output file
%  READ_ADCIRC_FORT62 reads in the contents of an ADCIRC station 
%  velocity file, typically output to IO unit 62 (fort.62).  
%  This output file is described in the ADCIRC Users Manual as 
%  "Depth-averaged Velocity Time Series at Specified Velocity 
%  Recording Stations (fort.62)". Time is converted from ADCIRC time 
%  to Gregorian time.
%  Call as: D=read_adcirc_fort62(fname,asd);
%
% READ_ADCIRC_FORT63 read ADCIRC global elevation output file
%  READ_ADCIRC_FORT63 reads in the contents of an ADCIRC global 
%  elevation file, typically output to IO unit 63 (fort.63).  
%  This output file is described in the ADCIRC Users Manual as 
%  "Elevation Time Series at All Nodes in the Model Grid (fort.63)". 
%  Time is converted from ADCIRC time to Gregorian time.
%  Call as: D=read_adcirc_fort63(fname,asd);
%
% READ_ADCIRC_FORT64 read ADCIRC global velocity output file
%  READ_ADCIRC_FORT64 reads in the contents of an ADCIRC global 
%  velocity file, typically output to IO unit 64 (fort.64).  
%  This output file is described in the ADCIRC Users Manual as 
%  "Depth-averaged Velocity Time Series at  All Nodes in the Model Grid  
%  (fort.64)". Time is converted from ADCIRC time 
%  to Gregorian time.
%  Call as: D=read_adcirc_fort64(fname,asd);
