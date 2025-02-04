%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% blochODE45
% This uses Matlab's ODE solver ode45, using the same interface as the RK4
% implementation. It is very slow because it calculates 50 intermediate
% points for each time step. I'm sure there's a faster way, but this
% implementation is really just to validate and compare the RK4 code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the time evolution of magnetization using the Bloch equations.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mnext = blochODE45(Minit, B1x, B1y, Omega, R1, R2, deltaT)
% The ode45 code requires that the first two parameters are time and the
% vector M (they use 'y' by default). However, we don't need t.
[t, Msim] = ode45(@(t,M) bloch_eqn(t, M, B1x, B1y, Omega, R1, R2), [0 deltaT], Minit.');

% It simulated many points, only need the last one
Mnext = Msim(end,:);
return;



function dMdt = bloch_eqn(t, M, b1x, b1y, Omega, R1, R2)
dMdt = [0 0 0]; % Initize the 3D vector

% The three components of dM/dt: 
dMdt(1) = (Omega * M(2)) - (b1y * M(3)) - (M(1) * R2);
dMdt(2) = (b1x * M(3)) - (Omega * M(1)) - (M(2) * R2);
dMdt(3) = (b1y * M(1)) - (b1x * M(2)) + ((1.0 - M(3)) * R1);

dMdt = dMdt';
return;



