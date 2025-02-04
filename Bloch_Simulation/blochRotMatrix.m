%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% blochRotMatrix
% Calculates the time evolution of magnetization described by the Bloch
% equations, using the rotation matrix approximation. 
% 
% This is really just for demonstration purposes - the RK4 variant is more
% accurate and should be preferred.
%
% Inputs:
%  	B1x,B1y - The x and y components of the complex RF amplitude in the
%       rotating frame. A purely amplitude modulated pulse will have all B1Y=0.
%		Units are rad/s [real].
%	deltaT - time step, in seconds [real]
% 	Omega - The off-resonance frequency, in rad/s/ [real]
%	R1 - Spin-lattice relaxation rate (1/T1), seconds^-1 [real]
%	R2 - Spin-spin relaxation rate (1/T2), seconds^-1 [real]
%	Minit - initial magnetization, normalized to 1 (ie equilibirum
%       magnetization M0=1) [Row vector of 3 reals]
%
% Outputs:
%	Mnext - The resultant magnetization. [Row vector of 3 reals]
%
% To calculate this over an RF pulse or pulse sequence, call it iteratively
% for each timestep and for each off-resonance component
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2020 Patrick Bolan, University of Minnesota
% I release this code into the public domain, without conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mnext = blochRotMatrix(Minit, B1x, B1y, Omega, R1, R2, deltaT)

% First, apply the RF. Convert to polar
B1amp = sqrt(B1x^2 + B1y^2);
B1phi = atan2(B1y, B1x);
alpha = -B1amp * deltaT;

% Apply this as Rz(phi) Rz(alpha) Rz(-phi)
Rzneg = [cos(B1phi), sin(B1phi), 0; -sin(B1phi), cos(B1phi), 0; 0, 0, 1];
Rzpos = [cos(B1phi), -sin(B1phi), 0; sin(B1phi), cos(B1phi), 0; 0, 0, 1];
Rx = [1, 0, 0; 0, cos(alpha), -sin(alpha); 0, sin(alpha), cos(alpha)];

Mrot = Rzpos * (Rx * (Rzneg * Minit.'));

% Then apply the precession and relaxation
OT = -Omega * deltaT;
Rprec = [cos(OT), -sin(OT), 0; sin(OT), cos(OT), 0; 0, 0, 1];
E1 = exp(-deltaT * R1);
E2 = exp(-deltaT * R2);

Mnext = [E2, 0, 0; 0, E2, 0; 0, 0, E1] * (Rprec * Mrot) + [0; 0; 1-E1];



