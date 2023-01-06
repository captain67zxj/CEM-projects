close all; clear; clc;
%%
a           =   1.6; % Radius of circle in wavelength
N           =   15*round(2*pi*a); % Number of segments
phi_i       =   180; % Angle of incidence in degree
%%
[Z_TM,I_TM,RCS_TM,Z_TE,I_TE,RCS_TE,phi]=RCS(a,N,phi_i);
%%